#!/bin/zsh

export PATH=$PATH:~/iInstalled/sratoolkit.3.0.5-mac64/bin

# set pwd to ERR1461 subdir
while read SRR  ; do
        echo $SRR
        numHere=$(ls -1 | grep $SRR | wc -l)
        if [[ $numHere -gt 0 ]]
        then
                echo "Already have a file for id: $SRR"
        else
                echo "Download $SRR"
                fasterq-dump $SRR -v --split-3 --outdir . 
                gzip ${SRR}.fastq
        fi
# done < ../subsetERR1459_SRR_Acc_List.txt
# done < ../subsetERR1460_SRR_Acc_List.txt
done < ../subsetERR1461_SRR_Acc_List.txt
# done < ../subsetERR1462_SRR_Acc_List.txt



# Round1
# while read SRR  ; do
#         echo $SRR
#         fasterq-dump $SRR -v --split-3 --outdir . 
# #done < SRR_Acc_List.txt

# Round1
# while read SRR  ; do
#         echo $SRR
#         fasterq-dump $SRR -v --split-3 --outdir . 
# #done < SRR_Acc_List_round2.txt

# Having >2k files in one folder made it hard to see the scripts in the viewer,
# not sure if it might have caused any other issues. (?)
# So I split the files into subdirs.
# I found that I had some but not all of files I should have in each subdir.
# So I split the Acc_List into subsets, so I can go through each subset in each folder.


# mkdir ERR1459
# mv ERR1459*fastq* ERR1459 
# cat SRR_Acc_List.txt | grep ERR1459 > subsetERR1459_SRR_Acc_List.txt

#wc -l subsetERR14*               
#>     679 subsetERR1459_SRR_Acc_List.txt
#>    1000 subsetERR1460_SRR_Acc_List.txt
#>    1000 subsetERR1461_SRR_Acc_List.txt
#>     400 subsetERR1462_SRR_Acc_List.txt
#>    3079 total

# ls -1 ERR1459 | wc -l
#      541
# ls -1 ERR1460 | wc -l                      
#      790
# ls -1 ERR1461 | wc -l                      
#      783
# ls -1 ERR1462 | wc -l                      
#      297


# First time in the ERR1459 subdir, I see this error:
# Failed to call external services.

# from https://hpc.nih.gov/apps/sratoolkit.html
# If you are getting this error while using fastq-dump and fasterq-dump: 
# Failed to call external services 
# try this: % mv ~/.ncbi ~/.ncbi.OLD
#
# I treid that. It did not help.
# Then I ran updates on my computer, restarted, turned my hotspot off and back on.... all the off and on again.
# Now its good. :)

