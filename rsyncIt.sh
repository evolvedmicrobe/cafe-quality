#!/bin/sh
rsync --verbose  --progress --stats --compress --rsh=/usr/bin/ssh \
      --recursive --times --perms --links --delete \
      --exclude "*.exe" --exclude "*.dll" --exclude "*.mdb" --exclude "*.csv" \
      --exclude "*.pdf" --exclude "*fasta.gz" --exclude '*pptx' --exclude 'data' \
      --exclude 'lib' --exclude '*.o' --exclude '*.d' --exclude '*.a' --exclude '*.so' \
      --exclude '*.dylib' --exclude '*bin/*' --exclude '*obj/*' \
      /Users/nigel/git/cafe-quality/* ndelaney@mp-f052:/home/UNIXHOME/ndelaney/git/cafe-quality