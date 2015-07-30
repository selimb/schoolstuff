@echo off
::http://superuser.com/questions/246837/how-do-i-add-text-to-the-beginning-of-a-file-in-bash
::START FROM mech_530
::Use %1 to refer to command line argument (the html)
::Use %~n1 to strip extension
set str=composites_assignment_6
echo %str%
set assnum=6
echo %assnum%
set notebookone=composites_assignment_6-design_1
set notebooktwo=composites_assignment_6-design_2
set notebookthree=composites_assignment_6-design_3
ipython nbconvert --to html --template html_nocode.tpl %notebookone%.ipynb
ipython nbconvert --to html --template html_nocode.tpl %notebooktwo%.ipynb
::ipython nbconvert --to html --template html_nocode.tpl %notebookthree%.ipynb
set notebookpdfone=notebook_%assnum%_one
set notebookpdftwo=notebook_%assnum%_two
set notebookpdfthree=notebook_%assnum%_three
wkhtmltopdf --image-quality 300 %notebookone%.html docs/%notebookpdfone%.pdf
wkhtmltopdf --image-quality 300 %notebooktwo%.html docs/%notebookpdftwo%.pdf
wkhtmltopdf --image-quality 300 %notebookthree%.html docs/%notebookpdfthree%.pdf
del %notebookone%.html
del %notebooktwo%.html
del %notebookthree%.html
cd docs
pdftk %notebookpdfone%.pdf %notebookpdftwo%.pdf %notebookpdfthree%.pdf cat output notebooks_6.pdf
del %notebookpdfone%.pdf
del %notebookpdftwo%.pdf
del %notebookpdfthree%.pdf
:: set /P thedate=Enter due date without year.
set thedate=November 18
set titlefile=titlepage_%assnum%
::PAUSE
echo \def\assnum{%assnum%}\def\thedate{%thedate%, 2014} | cat - composites_template.tex > %titlefile%.tex
::PAUSE
pdflatex %titlefile%.tex
pdftk %titlefile%.pdf notebooks_6.pdf ass6\ass6main.pdf cat output composites_%assnum%.pdf
del %titlefile%*.*