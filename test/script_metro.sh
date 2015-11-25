# script_metro.sh
# bash script to run multiple jobs on Sloan servers
# Authored by Arthur Delarue on 10/28/15

for fileName in `cat JSONparams/paramList.txt`
do
	JSONfileName="JSONparams/"$fileName
	outputFileName="Log_files/out_"${fileName%.json}".log"
	if [ -f $outputFileName ];
	then
		rm $outputFileName
	fi
	qsub -pe mthread 2 -l mem_free=4G -o $outputFileName ./juliaCallMetro.sh $JSONfileName
done

