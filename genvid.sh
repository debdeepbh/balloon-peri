#!/bin/bash
RNDN=$(date +%F-%H-%M-%S)-$(cat /dev/urandom | tr -cd 'a-z0-9' | head -c 5)
#ffmpeg -i img/tc_%5d.png -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2,zoompan=d=2" -r 25 -crf 30 img/${RNDN}.mp4
ffmpeg -i img/tc_%5d.png -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -r 25 -crf 30 img/${RNDN}.mp4
