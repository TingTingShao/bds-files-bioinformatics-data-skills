# Setting up and Managing a Bioinformatics Project

<hr>

## **NOTE FROM BOOK**<br />

<hr>

## project Directories and Directory structures

    >mkdir -p zmays-snps/{data/seqs,script,analysis}      
- zmays-snps: project name<br>
- data: all raw and intermediate data<br>
- data/seqs: dta processing as separate subjects<br>
- script: contain project wide scripts, very often divided into different moduler scripts<br>
- analysis: keep high level analyses, such as figures, and tables<br>
**each folder contains a `README.md` file**: what is in the repository and how it got there<br>

## PROJECT DOCUMENTATION
- document method and workflow: any command produces results, as well default values<br>
- document the origin of all data in the project directory: where, when, what has done to<br>
- document when you downloaded the data<br>
- record data version information<br>
- describe how you download the data<br>
- document the versions of the software you run, including the packages<br>

## USE DIRECTORIES TO DIVE UP YOUR PROJECT INTO SUBPROJECTS
- create directories to logically separate subprojects(for example data quality improvement, aligning, analyzing alignment results, etc.)<br>
## ORGANIZING DATA TO AUTOMATE FILE PROCESSING TASKS
    >touch seqs/zmays{A,B,C}_R{1,2}.fastq
**Numbering the files as 001, 002...014**<br>

## MARKDOWN FOR PROJECT NOTEBOOKS

## Markdown Editors

If you prefer to edit Markdown files directly in an WYSIWYG editor, here are
some good editors you might want to try:

 - [Mou](http://25.io/mou/) (OS X)
 - [MarkdownPad](http://markdownpad.com/) (Windows)

Classic editors like Vim and Emacs also have syntax highlighting for Markdown.

## Converting Markdown to PDF using Pandoc and LaTeX

Pandoc can also convert Markdown to PDF through a tool called LaTeX. If you
haven't heard of LaTeX, it's a powerful document markup language used in
mathematics and the sciences. I use LaTeX everyday in my work; it's the format
I (and many others) use for manuscripts. If you'd like to learn LaTeX, I
suggest the [LaTeX WikiBook](http://en.wikibooks.org/wiki/LaTeX). If you're
just looking to render HTML to PDF with Pandoc, you'll need to install LaTeX
first.

If you're on a Linux machine, you can install a LaTeX distribution through:

    $ sudo apt-get install texlive-full

If you're on a Mac, you can install a LaTeX distribution through
[MacTeX](https://tug.org/mactex/) (beware though; this is an exceptionally
large file).

Then, use:

    $ pandoc --from markdown --to latex notebook.md --output output.pdf

Note that we need to use `--output` rather than redirection, as this helps
Pandoc figure out you want a PDF file, not the intermediate LaTeX file.

## More on Wildcards

To prevent confusion, I've only taught the most common shell wildcards in this
chapter. There are two others you may be interested in. First, in Bash (this
does not working other shells like zsh), +[!BC]+ matches any character _other_
than B or C. Second, because the digits 0-9 are in order in ASCII, you can also
use character class ranges to specify single digit numbers like
+files-[3-7].txt+. However, remember the warning in this chapter that this does
_not_ work with double digit numbers like +files-10-14.txt+.

## Argument List Too Long

We'll discuss the solution to this in Chapter 12 (when we cover +xargs+), but a
quick note about the technical I mentioned in this chapter: this is not a
limitation of your shell, but actually a lower system call called +exec()+. If
these nerdy details interest you, here's [a good
explanation](http://stackoverflow.com/questions/4185017/maximum-number-of-bash-arguments-max-num-cp-arguments).

