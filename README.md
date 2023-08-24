HST1PASS
--------

Copied over from Jay's host website:

(https://www.stsci.edu/~jayander/HST1PASS/)[https://www.stsci.edu/~jayander/HST1PASS/]

on 2023-08-24 and posted to GitHub for easier install. To fully install the package:

1. git clone the repository
2. use the following two statements to fetch the LIB and USE_CASES directories. Warning: These are large.

wget -e robots=off --no-parent --no-host-directories -R "index.html*" -l 0 -r https://www.stsci.edu/~jayander/HST1PASS/LIB/

wget -e robots=off --no-parent --no-host-directories -R "index.html*" -l 0 -r https://www.stsci.edu/~jayander/HST1PASS/USE_CASES/

3. Compile the code as needed.

