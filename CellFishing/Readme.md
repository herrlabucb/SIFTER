# CellFishing README


The standalone source code and the real datasets are all included.
To use the code, software R is required and can be downloaded from: https://www.r-project.org\

It has been tested on a "normal" desktop, with versions listed below:
      platform       x86_64-apple-darwin15.6.0   
      arch           x86_64                     
      os             darwin15.6.0               
      system         x86_64, darwin15.6.0       
      status                                    
      major          3                          
      minor          6.1                        
      year           2019                       
      month          07                          
      day            05                          
      svn rev        76782                       
      language       R                           
      version.string R version 3.6.1 (2019-07-05)
      nickname       Action of the Toes

No other non-standard hardware or software is needed. Some R libraries is required and can be easily included as outlined in the source code.\
There is no need to install the code. The algorithm has been reported previously in a published research paper: 
Liu, K. et al. GeneFishing to reconstruct context specific portraits of biological processes. Proc. Natl. Acad. Sci. U. S. A. 116, 18943-18950 (2019).

To run the code with the data, you only need to run the one file: Fishing.R and the expected output is the file of format .RData, which consists of all the Capture Frequency for each cell. With the attached dataset it is expected to be finished in around 15 minutes.
