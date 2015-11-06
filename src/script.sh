# script.sh
# bash script to run multiple jobs on Sloan servers
# Authored by Arthur Delarue on 10/28/15

for fileName in `cat JSONparams/paramList.txt`
do
	JSONfileName="JSONparams/"$fileName
	outputFileName="out_"${fileName%.json}".log"
	if [ -f $FILE ];
	then
		rm $outputFileName
	fi
	qsub -pe mthread 10 -l mem_free=40G -o $outputFileName ./juliaCall.sh $JSONfileName
done