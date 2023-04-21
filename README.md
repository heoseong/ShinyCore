# ShinyCore
ShinyCore is available at https://stevenkimcsumb.shinyapps.io/ShinyCore, and it takes five steps. A file should be uploaded in Step 1. Wait for the message “Upload complete.” A csv file is recommended. The separator should be specified in Step 2. The default setting is the comma separator (,). In the data, the samples (accessions) should be organized by column or by row. It should be specified in Step 3. The default setting is by Column. See Table 1 and 2 for the two acceptable data structures. The maximum coverage (%) and the minimum coverage (%) should be specified in Step 4 and 5, respectively. The default values are 99 and 98, respectively. Then click on the Submit button. 

Sometimes, the server is not stable, or the server does not support a big dataset. In such a case, it is recommended to download the zipped folder <a href="https://github.com/heoseong/ShinyCore/blob/main/ShinyCore%202023-04-20.zip">at this link</a>. The user’s computer must have RStudio installed. Unzip the folder, and click on ShinyCore.R inside the folder. If the shiny package has not been installed on the user’s computer, install the shiny package. Click on the Run App button. The applet will be opened. The user’s guidelines are the same as given in the above paragraph.

The two wheat datasets are available in the subfolder named Data. As the entire collection is bigger, it takes a longer time to see the first entry. For example, using a Core i5 computer (the 10th Generation, 16.0 GB), it takes about 0.11 minutes (7 seconds) to see the first entry for the Wheat 1 dataset and 0.18 minutes (11 seconds) for the Wheat 2 dataset. 

Table 1. Samples (accessions) by column


| |Sample 1|Sample 2|...|Sample *n*|
|---|---|---|---|---|
|marker 1|1|2| | |
|marker 2|0|1| | |
|...| | | | |
|marker *m*|NA|0| | |

Table 2. Samples (accessions) by row


| |Marker 1|Marker 2|...|Marker *m*|
|---|---|---|---|---|
|Sample 1|G|T| | |
|Sample 2|A|C| | |
|...| | | | |
|Sample *n*|N|C| | |
