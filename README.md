# SVRecall
Used to help verify sample-unique SV calls after performing merging.  It reruns Sniffles on each of the panel samples focused on just the reads around the variant site, and if anything is close to the variant site (based on Jasmine distance), it is flagged and an updated entry is output.  


Note: Still very much a work in progress!
