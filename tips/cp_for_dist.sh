#!/bin/bash

if [ "$#" -eq 0 ]; then
  echo "エラー: 引数が必要です。" >&2
  exit 1
fi

if [ ! -d "./tips" ]; then
  echo "エラー: 'tips' ディレクトリが存在しません。" >&2
  exit 1
fi

mkdir $1
cp -r ../abmptools $1/
rm -rf $1/abmptools/build
rm -rf $1/abmptools/dist
rm -rf $1/abmptools/ABMPTools.egg-info
rm -rf $1/abmptools/abmptools/f90/src
rm -rf $1/abmptools/Makefile
rm -rf $1/abmptools/.git
