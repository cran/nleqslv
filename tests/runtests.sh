
for f in brdban chquad trig dslnex dslnexauto dslnexscaled
do  
	echo Running $f.R
	R CMD BATCH --arch=i386 --no-timing --no-save $f.R $f.Rout
done
	