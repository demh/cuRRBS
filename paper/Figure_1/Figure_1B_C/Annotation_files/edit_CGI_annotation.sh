#!/bin/bash

cp CGI.gtf CGI_mod.gtf

# Delete unused contigs

while read chrd; do
        echo $chrd
        sed -i.bak "/^$chrd/d" CGI_mod.gtf
done <chr_to_delete_CGI.txt

rm CGI_mod.gtf.bak

