#!/bin/bash
set -e

if `ls TAPE*_ex > /dev/null 2>&1`
then
  rm TAPE*_ex
fi

echo "\nRunning exec: lblrtm\n\n"

time ./lblrtm

if [ -e TAPE10 ] 
then
  mv TAPE10 TAPE10_ex
fi

if [ -e TAPE12 ] 
then
mv TAPE12 TAPE12_ex
fi

if [ -e TAPE27 ] 
then
mv TAPE27 TAPE27_ex
fi

if [ -e TAPE28 ] 
then
mv TAPE28 TAPE28_ex
fi

if [ -e TAPE6 ] 
then
mv TAPE6 TAPE6_ex
fi

if [ -e TMP* ] 
then
rm TMP*
fi
