library(shiny)

##########
### ui ###
##########

ui = fluidPage(

  titlePanel( "ShinyCore (Simple Version)" ),

  sidebarLayout(
    
    sidebarPanel(
      
      ### upload data
      
      fileInput( "target_upload", 
                 label = "Step 1: Upload the file",
                 accept = c( "text/csv", "text/comma-separated-values", ".csv" ) ),
      
      ### choose separator
      
      radioButtons( "separator", 
                    label = "Step 2: Choose the separator", 
                    choices = c( ",", ";", ":" ), 
                    selected = ",", inline = TRUE ),
      
      ### display data
      
      DT::dataTableOutput( "sample_table" ),
      
      ### data structure (sample by column or by row)
      
      selectInput( inputId = "sample.by", 
                   label = "Step 3: The samples are by", 
                   choices = c( "Column" = "column", "Row" = "row" ) ),
      
      ### input coverage
      
      numericInput( inputId = "coverage",
                    label = "Step 4: Enter the maximum coverage %",
                    value = 99,
                    min = 0,
                    max = 100 ),
      
      ### input coverage
      
      numericInput( inputId = "coverage.min",
                    label = "Step 5: Enter the minimum coverage %",
                    value = 98,
                    min = 0,
                    max = 100 ),
      
      ### submit button
      
      actionButton( "do", "Submit" )
      
    ),
    
    mainPanel(
      
      tabsetPanel(
        
        tabPanel( "Covering in Progress", verbatimTextOutput( "loc" ) ),
        tabPanel( "Coverage Plot", plotOutput( "coverage.plot", width = "450px", height = "450px" ) ),
        tabPanel( "Core Collection", verbatimTextOutput( "report1" ) ),
        tabPanel( "Core Evaluation", verbatimTextOutput( "report2" ) ),
        tabPanel( "Alive Counter", textOutput( "keepAlive" ) )
        
      )
      
    )
    
  ),

  #mainPanel( img( src = "logo.jpg", height = 200, width = 200 ) ),
  mainPanel( tags$a( "Instruction is available at https://github.com/heoseong/ShinyCore.", href = "https://github.com/heoseong/ShinyCore" ) ),
  mainPanel( h1("") ),
  
  ### keep alive!
  ### source: https://tickets.dominodatalab.com/hc/en-us/articles/360015932932-Increasing-the-timeout-for-Shiny-Server
  
  tags$head(
    HTML(
      "
          <script>
          var socket_timeout_interval
          var n = 0
          $(document).on('shiny:connected', function(event) {
          socket_timeout_interval = setInterval(function(){
          Shiny.onInputChange('count', n++)
          }, 15000)
          });
          $(document).on('shiny:disconnected', function(event) {
          clearInterval(socket_timeout_interval)
          });
          </script>
          " ) ),
  
  ### suppress error message!
  
  tags$style( type = "text/css",
              ".shiny-output-error { visibility: hidden; }",
              ".shiny-output-error:before { visibility: hidden; }" ) )

##############
### server ###
##############

server = function(input, output, session ) {
  
  ### data size: up to 30 MB
  
  options( shiny.maxRequestSize = 30 * 1024 ^ 2 )
  
  g = reactiveValues()
  
  ### execute (when action button is pressed)
  
  observeEvent( input$do, {
    
    start = proc.time()
    
    f.count = function(x) { length( unique( x[ !is.na(x) ] ) ) }
    
    f.count.na = function(x) { sum( is.na(x) ) }
    
    f.numeric = function(x) { as.numeric( as.factor(x) ) }
    
    f.shannon = function(x) { 
      
      p = as.numeric( prop.table( table(x) ) )
      - sum( log( p ^ p, max( length(p), 2 ) ) )
      
    }
    
    #############
    ### READY ###
    #############
    
    ### read in data
    
    df_products_upload = reactive( {
      
      inFile = input$target_upload
      if (is.null(inFile))
        return(NULL)
      df = read.csv( inFile$datapath, header = TRUE, sep = input$separator )
      return(df)
      
    } )
    
    ### display the data
    
    #output$sample_table = DT::renderDataTable({
    
    #df = df_products_upload()
    #DT::datatable(df)
    
    #})
    
    ### format the data
    
    data0 = df_products_upload()
    
    if ( input$sample.by == "column" ) {
      
      c.names = colnames( data0[,-1] ) ### column names
      r.names = data0[,1] ### row names
      data = t( data0[,-1] ) ### transpose the data after removing the first column
      data = data.frame(data) ### data frame
      rownames(data) = c.names
      colnames(data) = r.names
      count.temp = as.numeric( unlist( lapply( data, f.count ) ) )
      data = data[ , count.temp >= 2 ] ### include markers with two or more genotypes
      
    }
    
    if ( input$sample.by == "row" ) {
      
      c.names = colnames( data0[,-1] ) ### column names
      r.names = data0[,1] ### row names
      data = data0[,-1]
      data = data.frame(data) ### data frame
      rownames(data) = r.names
      colnames(data) = c.names
      count.temp = as.numeric( unlist( lapply( data, f.count ) ) )
      data = data[ , count.temp >= 2 ] ### include markers with two or more genotypes
      
    }
    
    ### order markers by rarity
    
    temp = lapply( data, table )
    temp = lapply( temp, prop.table )
    temp = as.numeric( unlist( lapply( temp, min ) ) )
    data = data[ , order( temp, decreasing = FALSE ) ]
    
    ### data dimension
    
    n = dim(data)[1] ### number of rows (samples)
    k = dim(data)[2] ### number of columns (markers)
    sample.name = rownames(data)
    marker.name = colnames(data)
    
    ### count genotypes per marker
    
    count.all = as.numeric( sapply( data, f.count ) )
    
    ###############
    ### ROUND 1 ###
    ###############
    
    data.noncore = data
    
    ### choose the sample with the minimum number of NA's
    
    temp = as.numeric( apply( data.noncore, 1, f.count.na ) )
    i.core = which.min(temp) ### initial core set
    
    data.core.new = data[i.core,]
    data.core = data.core.new
    data.noncore = data.noncore[-i.core,]
    count.core = as.numeric( sapply( data.core, f.count ) )
    marker.covered = which( count.core == count.all )
    
    ### report
    
    sample.new = rownames(data.core.new)
    k0 = length(count.all)
    coverage.new = mean( count.core / count.all ) * k0 / k + ( k - k0 ) / k 

    now = proc.time()
    time.elapsed = ( now - start )[3] / 60
    line0 = paste0( "iteration = ", 1, 
                    "; time = ", formatC( time.elapsed, format = "f", digit = 3 ), 
                    " min; CV = ", formatC( coverage.new * 100, format = "f", digit = 1 ), 
                    "%; new sample = ", sample.new )
    
    ### remove markers that are covered
    
    if ( length(marker.covered) > 0 ) { 
      
      data.core = data.core[,-marker.covered]
      data.noncore = data.noncore[,-marker.covered]
      data.core.new = data.core.new[,-marker.covered]
      count.all = count.all[-marker.covered]
      count.core = count.core[-marker.covered]
      
    }
    
    ### assign NA to types that are covered
    
    data.noncore[ data.noncore == sapply( data.core.new, rep, dim(data.noncore)[1] ) ] = NA
    
    ###############
    ### ROUND 2 ###
    ###############
    
    ### start here!
    
    rv = reactiveValues( lines = character(0), i = 1 )
    
    ### carry important results to reactive values
    
    rv$data.core = data.core
    rv$data.noncore = data.noncore
    rv$data.core.new = data.core.new
    
    rv$count.all = count.all
    rv$count.core = count.core
    
    rv$coverage.new = coverage.new
    rv$coverage.history = coverage.new
    rv$sample.history = sample.new
    
    rv$lines = line0
    
    rv$time.elapsed = 0
    rv$pause = 0
    
    observe( {
      
      # Re-execute this reactive expression immediately after 
      # it finishes (though you can also specify a longer time:
      # e.g. 1000 means re-execute after 1000 milliseconds have
      # passed since it was last executed)
      
      invalidateLater( 0, session )
      
      isolate( {
        
        ### thickening
        
        if ( rv$coverage.new >= input$coverage / 100 & rv$pause == 0 ) { 
          
          now = proc.time()
          time.elapsed = ( now - start )[3] / 60
          T1 = formatC( time.elapsed, format = "f", digits = 3 )
          
          core1 = rv$sample.history
          N1 = length(core1)
          N0 = min( which( rv$coverage.history > input$coverage.min / 100 ) )
          
          odds = function(x) { p = prop.table( table(x) ); ( 1 - p ) / p }
          
          score.sheet = data
          k = dim(data)[2]
          
          for (j in 1:k ) { score.sheet[,j] = as.factor( score.sheet[,j] ) }
          temp = sapply( score.sheet, odds )
          
          for (j in 1:k ) {
            
            temp.name = names( temp[[j]] )
            temp.score = as.numeric( temp[[j]] )
            score.sheet[,j] = temp.score[ match( score.sheet[,j], temp.name ) ]
            
          }
          
          score.sample = apply( score.sheet, 1, sum, na.rm = TRUE )
          score.sample = score.sample[ order( score.sample, decreasing = TRUE ) ]
          order.sample = names(score.sample)[ !names(score.sample) %in% core1[1:N0] ]
          core2 = c( core1[1:N0], order.sample[ 1:( N1 - N0 ) ] )
          
          data.temp1 = data[ rownames(data) %in% core1, ]
          data.temp2 = data[ rownames(data) %in% core2, ]
          
          ### RS
          
          RS1 = log10( sum( score.sample[ names(score.sample) %in% core1 ] ) )
          RS1 = formatC( RS1, format = "f", digits = 3 )
          
          RS2 = log10( sum( score.sample[ names(score.sample) %in% core2 ] ) )
          RS2 = formatC( RS2, format = "f", digits = 3 )
          
          ### SH
          
          SH1.markers = as.numeric( apply( data.temp1, 2, f.shannon ) )
          SH1 = mean(SH1.markers)
          SH1 = formatC( SH1, format = "f", digits = 3 )
          
          SH2.markers = as.numeric( apply( data.temp2, 2, f.shannon ) )
          SH2 = mean(SH2.markers)
          SH2 = formatC( SH2, format = "f", digits = 3 )
          
          ### CV
          
          CV1 = mean( as.numeric( sapply( data.temp1, f.count ) ) / count.all )
          CV1 = formatC( CV1, format = "f", digits = 3 )
          
          CV2 = mean( as.numeric( sapply( data.temp2, f.count ) ) / count.all )
          CV2 = formatC( CV2, format = "f", digits = 3 )
          
          ### MR
          
          n.temp = length(core1)
          MR1.box = matrix( NA, nrow = n.temp, ncol = n.temp )
          diag(MR1.box) = NA
          colnames(data.temp1) = NULL; rownames(data.temp1) = NULL
          data.temp1 = as.matrix(data.temp1)

          for (i in 2:n.temp ) {
            
            if ( i == 2 ) mr.temp = sqrt( mean( c( data.temp1[1,] ) != c( data.temp1[i,] ), na.rm = TRUE ) )
            if ( i >= 3 ) mr.temp = sqrt( apply( t( data.temp1[1:(i - 1),] ) != c( data.temp1[i,] ), 2, mean, na.rm = TRUE ) )
            MR1.box[1:(i - 1),i] = mr.temp

          }
          
          MR1 = mean( MR1.box, na.rm = TRUE )
          MR1 = formatC( MR1, format = "f", digits = 3 )
          
          n.temp = length(core2)
          MR2.box = matrix( NA, nrow = n.temp, ncol = n.temp )
          colnames(data.temp2) = NULL; rownames(data.temp2) = NULL
          data.temp2 = as.matrix(data.temp2)
          
          for (i in 2:n.temp ) {
            
            if ( i == 2 ) mr.temp = sqrt( mean( c( data.temp2[1,] ) != c( data.temp2[i,] ), na.rm = TRUE ) )
            if ( i >= 3 ) mr.temp = sqrt( apply( t( data.temp2[1:(i - 1),] ) != c( data.temp2[i,] ), 2, mean, na.rm = TRUE ) )
            MR2.box[1:(i - 1),i] = mr.temp
            
          }
          
          MR2 = mean( MR2.box, na.rm = TRUE )
          MR2 = formatC( MR2, format = "f", digits = 3 )
          
          now = proc.time()
          time.elapsed = ( now - start )[3] / 60
          T2 = formatC( time.elapsed, format = "f", digits = 3 )
          
          report1 = data.frame( Before.Thickening = core1, 
                                After.Thickening = core2 )
          
          report2 = data.frame( Before.Thickening = c( CV1, SH1, MR1, RS1, T1 ),
                                After.Thickening = c( CV2, SH2, MR2, RS2, T2 ) )
          
          rownames(report2) = c( "CV", "SH", "MR", "RS", "Time.Minutes" )
          
          output$report1 = renderPrint(report1)
          output$report2 = renderPrint(report2)
          
          rv$pause = 1
          
        } 
        
        ### no more updates
        
        if ( rv$coverage.new >= input$coverage / 100 & rv$pause == 1 ) { return() } 
        
        ### update core
        
        else {
          
          ### choose the sample with the minimum number of NA's
          
          rv$i = rv$i + 1
          temp = as.numeric( apply( rv$data.noncore, 1, f.count.na ) )
          i.core = which.min(temp)
          rv$data.core.new = rv$data.noncore[i.core,]
          
          rv$data.core = rbind( rv$data.core, rv$data.core.new )
          rv$data.noncore = rv$data.noncore[-i.core,]
          rv$count.core = as.numeric( sapply( rv$data.core, f.count ) )
          
          marker.covered = which( rv$count.core == rv$count.all )
          
          ### report
          
          rv$sample.new = sample.new = rownames(rv$data.core.new)
          k0 = length(rv$count.all)
          rv$coverage.new = mean( rv$count.core / rv$count.all ) * k0 / k + ( k - k0 ) / k 
          
          now = proc.time()
          rv$time.elapsed = time.elapsed = ( now - start )[3] / 60
          rv$lines = rbind( rv$lines, 
                            paste0( "iteration = ", rv$i, 
                                    "; time = ", formatC( time.elapsed, format = "f", digit = 3 ), 
                                    " min; CV = ", formatC( rv$coverage.new * 100, format = "f", digit = 1 ), 
                                    "%; new sample = ", sample.new ) )
          colnames(rv$lines) = ""
          
          coverage.history = rv$coverage.history = c( rv$coverage.history, rv$coverage.new )
          sample.history = rv$sample.history = c( rv$sample.history, rv$sample.new )

          ### remove markers that are covered
          
          if ( length(marker.covered) > 0 ) { 
            
            rv$data.core = rv$data.core[,-marker.covered]
            rv$data.noncore = rv$data.noncore[,-marker.covered]
            rv$data.core.new = rv$data.core.new[,-marker.covered]
            rv$count.all = rv$count.all[-marker.covered]
            rv$count.core = rv$count.core[-marker.covered]
            
          }
          
          ### assign NA to types that are covered
          
          rv$data.noncore[ rv$data.noncore == sapply( rv$data.core.new, rep, dim(rv$data.noncore)[1] ) ] = NA
          
        } ### end of else
        
      } ) ### end of isolate
      
    } ) ### end of observe
    
    ### print the updated core
    
    output$loc = renderPrint(rv$lines)
    
    ### figure for coverage growth per iteration
    
    output$coverage.plot = renderPlot( {
      
      n.core = length(rv$coverage.history)
      plot( 0:n.core, c( 0, rv$coverage.history ), type = "l", 
            col = 4, lty = 1, 
            ylim = c(0,1), xlim = c(1,n.core), axes = FALSE,
            xlab = "Iteration", ylab = "Coverage", main = "Covering Phase" )
      axis( 1, c( -100, 0:n.core ) )
      axis( 2, seq( -0.2, 1, 0.2 ) )
      abline( h = input$coverage / 100, col = 8, lty = 3 )

      legend( "bottomright", col = c(4,8), lty = c(1,3), 
              legend = c( "CV", "Goal" ), bty = "n", cex = 1 )
      
    } ) ### end of renderPlot
    
  } ) ### end of observeEvent
  
  ### keep alive!
  
  output$keepAlive = renderText( { req(input$count); paste( "keep alive ", input$count ) } )
  
}

#############
### shiny ###
#############

shinyApp( ui, server )