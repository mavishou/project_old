#!/bin/bash

function run_cuffquant {
	sample=$1
	cmd="cuffquant -p 10 -o ../cuff/${sample}/ /lustre/user/houm/projects/AnnoLnc/V4expCaculate.combined.gtf ${sample}.bam"
	echo $cmd
	$cmd
}

while read s
do
	now_date=`date +%y%m%d`
	mkdir -p ../cuff/$s
	(echo "Running cuffquant for $s"
	time run_cuffquant $s) 2>&1 | tee ../cuff/logs/${now_date}-${s}_cuffquant.log
done