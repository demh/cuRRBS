#!/bin/bash

cp CGI_mouse_mm10.gtf CGI_mouse_mm10_mod.gtf

# Delete unused contigs

while read chrd; do
        echo $chrd
        sed -i.bak "/^$chrd/d" CGI_mouse_mm10_mod.gtf
done <chr_to_delete.txt

rm CGI_mouse_mm10_mod.gtf.bak

