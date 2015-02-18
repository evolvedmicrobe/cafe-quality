#!/bin/bash

export bucket='Releases'
if [ "$1" = "--pre" ]; then
  export bucket='Candidates'
fi

export destdir=/cygdrive/s/Software/$bucket/Edna
export srcdir='.'
export build=`p4 changes -m 1 | cut -f 2 -d ' '`
export who=`whoami`



if [ "$build" == "${build%+}" ]; then
  echo $build
else
  echo "Warning: you have files checked out in your 'Default' changelist."
  echo "If these are not related to the build (e.g. Matlab or Mono or"
  echo "somesuch), please put them in a separate numbered changelist"
  echo "and make another build."
fi


echo build=$build
mkdir -p $destdir/$build
cp -r $srcdir/bin/Release/* $destdir/$build
mkshortcut.exe --name "Edna" --desc="Latest Edna ($build) built by $who" $destdir/$build/Edna.exe
cp "Edna.lnk" $destdir
