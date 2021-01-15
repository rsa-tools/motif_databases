# how this repo as been created

### Clone old and new empty repo

```
git clone git@depot.biologie.ens.fr:rsat ens_rsat_motif_databases
git clone git@github.com:rsa-tools/motif_databases.git github_motif_databases
```

### Use filter-repo to filter the directory

```
cd ens_rsat_motif_databases
git remote remove origin # no needed, but it is more secure (should not be done if we want to update this repo)
git filter-repo --subdirectory-filter public_html/motif_databases --force --preserve-commit-hashes
```
See https://github.com/newren/git-filter-repo


### Get the commit in the new repo
```
cd ..
cd github_motif_databases
git remote add ens_updated  ../ens_rsat_motif_databases
git fetch ens_updated
git merge ens_updated/master -s recursive -X theirs
```


### Use git lfs for large file

- config for tf file
```
git lfs install
git lfs track "*.tf"
git lfs migrate import --include="*.tf"
git add .gitattributes
git commit -am 'use git lfs for *.tf file'
```
See: https://git-lfs.github.com/

- import everything (for every branch and other git ref)
```
git lfs migrate info --everything # get info
git lfs migrate import --everything --include="*.tf"
```


### Push in new repo on main branch
```
git push -u -f origin master
```
