%%
%% This is file `ncsuthesis.cls',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% ncsuthesis.dtx  (with options: `class')
%% --------------------------------------------------------------------
%% ncsuthesis --- A thesis class conforming to NC State ETD guidelines.
%% E-mail: jmhykes@ncsu.edu
%% Released under the LaTeX Project Public License v1.3c or later
%% See http://www.latex-project.org/lppl.txt
%% --------------------------------------------------------------------
%% 
\usepackage{fancyhdr}
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{ncsuthesis}[2012/01/26 v0.4 NC State ETD conformant thesis class]
%% -------------------------------------------------------------------------- %%
%% ------------------------------ Thesis type ------------------------------- %%
%% -------------------------------------------------------------------------- %%
\newcommand{\ncsu@thesistype}{phd}
\DeclareOption{phd}{\renewcommand{\ncsu@thesistype}{phd}}
\DeclareOption{ms}{\renewcommand{\ncsu@thesistype}{ms}}

%% Spacing options
\newcommand{\ncsu@defaultspacing}{\onehalfspacing}
\DeclareOption{singlespacing}{
     \renewcommand{\ncsu@defaultspacing}{\singlespacing}}
\DeclareOption{onehalfspacing}{
     \renewcommand{\ncsu@defaultspacing}{\onehalfspacing}}
\DeclareOption{doublespacing}{
     \renewcommand{\ncsu@defaultspacing}{\doublespacing}}

\newenvironment{lscape}{\begin{lscapenum}}{\end{lscapenum}}
\DeclareOption{hardcopy}{
    \renewenvironment{lscape}{\begin{lscapenotnum}}{\end{lscapenotnum}}
}

\DeclareOption{draft}{\PassOptionsToClass{\CurrentOption}{book}}
\DeclareOption{10pt}{\PassOptionsToClass{\CurrentOption}{book}}
\DeclareOption{11pt}{\PassOptionsToClass{\CurrentOption}{book}}
\DeclareOption{12pt}{\PassOptionsToClass{\CurrentOption}{book}}

\DeclareOption*{\PackageWarning{ncsuthesis}{Unknown option `\CurrentOption'. Ignoring}}
\ExecuteOptions{phd,onehalfspacing} % the default option is <phd, onehalfspacing>
\ProcessOptions
\LoadClass[openany,oneside,titlepage,letterpaper]{book}

\RequirePackage{setspace}
\RequirePackage{hyphenat}
\RequirePackage{nolbreaks}
\RequirePackage{ifthen}
\RequirePackage{sectsty}
\RequirePackage{tocloft}
\RequirePackage{graphicx}
\RequirePackage{calc}

%%----------------------------------------------------------------------------%%
%%---------------------------- Formatting Options ----------------------------%%
%% -------------------------------------------------------------------------- %%
%% Capital letter height
\newlength{\ncsu@capitalheight}
\settoheight{\ncsu@capitalheight}{D}

%% Line height
\newlength{\ncsu@lineheight}
\setlength{\ncsu@lineheight}{2.0ex}

%% Leading height: The leading changes thoughout the document, so this
%%  length must be updated
\newlength{\ncsu@leading}

%% We must add two double spaces plus one capital-height (for the page number)
%%  to footskip.  Footskip is measured from the top of the bottom margin.
%%  (See geometry package documentation, note that includefoot option is used.)
\ncsu@defaultspacing
\newlength{\ncsu@footskip}
\setlength{\ncsu@footskip}{\ncsu@capitalheight}
\addtolength{\ncsu@footskip}{2\baselineskip}
\addtolength{\ncsu@footskip}{0.5ex}  % A little extra space to ensure there
                                     %  are 2 full double spaced lines
\RequirePackage[includefoot,
                margin=1in,
                lmargin=1.25in,
                footskip=\ncsu@footskip]{geometry}

%% -------------------------------------------------------------------------- %%
%% Disposition format -- any titles, headings, section titles
\newcommand{\ncsu@dispositionformat}{}
\allsectionsfont{\ncsu@dispositionformat}
\providecommand{\dispositionformat}[1]{
    \renewcommand{\ncsu@dispositionformat}{#1}
    \allsectionsfont{\ncsu@dispositionformat}
}

%% Heading/Chapter/Section format
\newcommand{\ncsu@headingformat}{\ncsu@dispositionformat}
\providecommand{\headingformat}[1]{
    \renewcommand{\ncsu@headingformat}{\ncsu@dispositionformat #1}
}

%% \chapter*-to-first line space
%% This applies to preliminary headings, and ToC, LoF, LoT
\newlength{\ncsu@afterschapsep}
\setlength{\ncsu@afterschapsep}{3\ncsu@lineheight}

%% -------------------------------------------------------------------------- %%
%% Headers/footers - page number
\pagestyle{fancy}
\fancyhf{}
\lhead{\nouppercase{\leftmark}}
\rhead{\thepage}

%% -------------------------------------------------------------------------- %%
%% Vertical space around tables and figures, triple spacing
\setlength{\textfloatsep}{3\ncsu@lineheight}
\setlength{\intextsep}{3\ncsu@lineheight}

%%----------------------------------------------------------------------------%%
%%---------------------------  Content Options -------------------------------%%
%% -------------------------------------------------------------------------- %%
%% Author name
\newcommand{\ncsu@authorfirstmiddle}{First Middle}
\newcommand{\ncsu@authorlast}{Last}

%% Command to test for period after middle name
\newboolean{@ncsu@needsperiod}
\newcommand{\ncsu@testlastchar}[2]{%
    \ifthenelse{\equal{#2}{\relax}}{
        \ifthenelse{\equal{#1}{.}}{%
            \setboolean{@ncsu@needsperiod}{false}
        }{% #1 != .
            \setboolean{@ncsu@needsperiod}{true}
        }
        \let\next=\relax
    }
    {% #2 != relax
        \let\next=\ncsu@testlastchar
    }
    \next#2
}

\providecommand{\student}[2]{
    \renewcommand{\ncsu@authorfirstmiddle}{#1}
    \renewcommand{\ncsu@authorlast}{#2}
    \ncsu@testlastchar#1\relax
    \author{\ncsu@authorfirstmiddle{} \ncsu@authorlast}
}

%% Degree year - in case the degree year doesn't match the current year.
\newcommand{\ncsu@year}{\the\year}
\providecommand{\degreeyear}[1]{
    \renewcommand{\ncsu@year}{#1}
}

%% -------------------------------------------------------------------------- %%
%% Thesis Title
\newcommand{\ncsu@thesistitle}{My Thesis Title}
\providecommand{\thesistitle}[1]{\renewcommand{\ncsu@thesistitle}{#1}}

%% -------------------------------------------------------------------------- %%
%% ---------------------  Chapter Commands ---------------------------------- %%
%% -------------------------------------------------------------------------- %%

%% -------------------------------------------------------------------------- %%
%%  Redefine the \chapter* driver \@makeschapterhead. This is used for the
%%   preliminary headings.
%%   The spacing between the title baseline and the top of capital letters on
%%   the first line is specified by \ncsu@afterschapsep. However, this vertical
%%   height built by two separate lengths:
%%       1. the space between the baseline of the title and the top of the
%%           first line (including leading above the line).
%%              (\ncsu@afterschapbeforeleadingsep)
%%       2. From the top of the first line leading to the top of
%%           capital letters on the first line. (\ncsu@leading)
\newlength{\ncsu@afterschapbeforeleadingsep}
\AtBeginDocument{% at begin doc is necessary since fncychap
                 %  redefines \@makeschapterhead
    \renewcommand{\@makeschapterhead}[1]{%
        % Measure the leading here (which depends on the current line spacing),
        % and adjust the spacing accordingly between the baseline of the
        % heading and the top of capital letters on the first line.
        \setlength{\ncsu@leading}{\baselineskip-\ncsu@capitalheight}%
        \setlength{\ncsu@afterschapbeforeleadingsep}%
                  {\ncsu@afterschapsep-\ncsu@leading}%
        {\noindent\centering%
          {\ncsu@headingformat #1}%
          \par\nobreak%
          \vspace{\ncsu@afterschapbeforeleadingsep}}%
    }
}

%% -------------------------------------------------------------------------- %%
%% Abstract page commands
%%  - hyphnat package \nohyphens prevents title from being hyphenated.
%%  - \nolbreaks removes any line breaks in title (but you must use \linebreak
%%     instead of \\
\newcommand{\abstractname}{Abstract}
\newenvironment{abstract}{%
\singlespacing
	\chapter*{\abstractname}
   \begin{center}
   \begin{large}
   		\textbf{\ncsu@thesistitle} \\[1.5em]
   \end{large}
   By \@author
   \end{center}
    
%   \ncsu@defaultspacing
   \thispagestyle{empty}
}%
{\endtitlepage \pagestyle{empty}}

%% -------------------------------------------------------------------------- %%
%% Abstract page commands
\newcommand{\makecopyrightpage}{
	\thispagestyle{empty}
    \begin{titlepage}
    \centering
      \quad
      \vspace{0.47\textheight} %SPB

      \copyright{} Copyright \ncsu@year{} by \@author \\[1em]
      All Rights Reserved

      \vfill
    \end{titlepage}
}

%% -------------------------------------------------------------------------- %%
%% Title page commands

\newcommand{\ncsu@committeesize}{4}
\providecommand{\committeesize}[1]{\renewcommand{\ncsu@committeesize}{#1}}

%% To allow for Chair vs Co-chairs:
\newcommand{\ncsu@chairtype}{}

%% For each committee member, we have their name and title.
%%  By default, their title is empty.  We use an optional argument
%%  to set the title of committee members or to set the chair's
%%  title to something other than `Chair of Advisory Committee'.
\newcommand{\ncsu@chair}{Advisor Q. Chair}
\newcommand{\ncsu@chairtitle}{Chair of Advisory Committee}
\providecommand{\chair}[2][]{
    \renewcommand{\ncsu@chair}{#2}
    % Checking if optional title argument is present.
    \ifthenelse{\equal{#1}{}}{}{\renewcommand{\ncsu@chairtitle}{#1}}
    % Declare that we have a chair:
    \renewcommand{\ncsu@chairtype}{chair}
}

\newcommand{\ncsu@cochairI}{Co-chair 1 of Committee}
\newcommand{\ncsu@cochairItitle}{Co-chair of Advisory Committee}
\providecommand{\cochairI}[2][]{
    \renewcommand{\ncsu@cochairI}{#2}
    \ifthenelse{\equal{#1}{}}{}{\renewcommand{\ncsu@cochairItitle}{#1}}
    % Declare that we have co-chairs.
    \renewcommand{\ncsu@chairtype}{cochair}
}

\newcommand{\ncsu@cochairII}{Co-chair 2 of Committee}
\newcommand{\ncsu@cochairIItitle}{Co-chair of Advisory Committee}
\providecommand{\cochairII}[2][]{
    \renewcommand{\ncsu@cochairII}{#2}
    \ifthenelse{\equal{#1}{}}{}{\renewcommand{\ncsu@cochairIItitle}{#1}}
    \renewcommand{\ncsu@chairtype}{cochair}
}

\newcommand{\ncsu@memberI}{Member 1 of Committee}
\newcommand{\ncsu@memberItitle}{}
\providecommand{\memberI}[2][]{
    \renewcommand{\ncsu@memberI}{#2}
    \ifthenelse{\equal{#1}{}}{}{\renewcommand{\ncsu@memberItitle}{#1}}
}

\newcommand{\ncsu@memberII}{Member 2 of Committee}
\newcommand{\ncsu@memberIItitle}{}
\providecommand{\memberII}[2][]{
    \renewcommand{\ncsu@memberII}{#2}
    \ifthenelse{\equal{#1}{}}{}{\renewcommand{\ncsu@memberIItitle}{#1}}
}

\newcommand{\ncsu@memberIII}{Member 3 of Committee}
\newcommand{\ncsu@memberIIItitle}{}
\providecommand{\memberIII}[2][]{
    \renewcommand{\ncsu@memberIII}{#2}
    \ifthenelse{\equal{#1}{}}{}{\renewcommand{\ncsu@memberIIItitle}{#1}}
}

\newcommand{\ncsu@memberIV}{Member 4 of Committee}
\newcommand{\ncsu@memberIVtitle}{}
\providecommand{\memberIV}[2][]{
    \renewcommand{\ncsu@memberIV}{#2}
    \ifthenelse{\equal{#1}{}}{}{\renewcommand{\ncsu@memberIVtitle}{#1}}
}

\newcommand{\ncsu@memberV}{Member 5 of Committee}
\newcommand{\ncsu@memberVtitle}{}
\providecommand{\memberV}[2][]{
    \renewcommand{\ncsu@memberV}{#2}
    \ifthenelse{\equal{#1}{}}{}{\renewcommand{\ncsu@memberVtitle}{#1}}
}

\newcommand{\ncsu@titlepageparagraph}{
 \ifthenelse{\equal{\ncsu@thesistype}{phd}}{
   A dissertation submitted to the Graduate Faculty of\\
   North Carolina State University\\
   in partial fulfillment of the\\
   requirements for the Degree of\\
   Doctor of Philosophy
 }{} % end ifthenels phd
 \ifthenelse{\equal{\ncsu@thesistype}{ms}}{
   A thesis submitted to the Graduate Faculty of\\
   North Carolina State University\\
   in partial fulfillment of the\\
   requirements for the Degree of\\
   Master of Science
 }{} % end ifthenels phd
}

\newcommand{\ncsu@program}{Degree Program}
\providecommand{\program}[1]{\renewcommand{\ncsu@program}{#1}}

\newcommand{\ncsu@citystate}{Raleigh, North Carolina}

\newcommand{\ncsu@signatureline}[1][]{
    \centering
     \rule{0.7\textwidth}{0.4pt}
     #1
}

\newcommand{\ncsu@committeesignatures}{
	\rule{0.7\textwidth}{1pt} \\
	\ncsu@memberI, Ph.D\\
	Advisor \\[2em] 
	\rule{0.7\textwidth}{1pt} \\
	\ncsu@memberII, Ph.D \\
	Committee Member \\[2em]
	\rule{0.7\textwidth}{1pt} \\
	\ncsu@memberIII, Ph.D \\
	Committee Member\\[2em]
	\rule{0.7\textwidth}{1pt} \\
	\ncsu@memberIV, Ph.D \\
	Committee Member\\[2em]
	\rule{0.7\textwidth}{1pt} \\
	\ncsu@memberV, Ph.D \\
	Committee Member\\[3em]
}

\newcommand{\makedistribution}{
\thispagestyle{empty}
\begin{LARGE}
\noindent\textbf{Distribution Agreement} \\[2em]
\end{LARGE}
In presenting this thesis or dissertation as a partial fulfillment of the requirements for an advanced degree from Emory University, I hereby grant to Emory University and its agents the non-exclusive license to archive, make accessible, and display my thesis or dissertation in whole or in part in all forms of media, now or hereafter known, including display on the world wide web.  I understand that I may select some access restrictions as part of the online submission of this thesis or dissertation.  I retain all ownership rights to the copyright of the thesis or dissertation.  I also retain the right to use in future works (such as articles or books) all or part of this thesis or dissertation. \\[3em]
Signature: \\[3em]
\begin{tabular}{lcl}
\rule{0.4\textwidth}{1pt} && \rule{0.2\textwidth}{1pt}\\
Xiang Cheng && Date
\end{tabular}
}

\newcommand{\makeapproval}{
     \thispagestyle{empty}
     \clearpage
     \singlespacing
      \begin{center}
        \ncsu@thesistitle \\[1.5em]
        by \\
        \@author \\
         Doctor of Philosophy \\
       \ncsu@program\\[3.0\ncsu@lineheight]
       \ncsu@committeesignatures
       
       Accepted: \\
       \vfill
       \vfill
       \rule{0.7\textwidth}{1pt} \\
		Lisa A. Tedesco, Ph.D.\\
		Dean of the James T. Laney School of Graduate Studies
		\vfill
		\vfill
       \rule{0.5\textwidth}{1pt} \\
       Date
      \end{center}
      \thispagestyle{empty}
     \ncsu@defaultspacing
}

\newcommand{\makecover}{
     \thispagestyle{empty}
     \clearpage
     \singlespacing
      \begin{center}
        \begin{Large}
        \textbf{\ncsu@thesistitle} \\[2em]
        \end{Large}
               
        by \\ [2em]
        \@author \\
        M.S., Marquette University, 2011 \\
        B.S., Sichuan University, 2009 \\[2em]
        Advisor: Stefan Boettcher, Ph.D. \\[5em]
        An abstract of \\
        A dissertation submitted to the Faculty of the \\
        James T. Laney School of Graduate Studies of Emory University \\
        in partial fulfillment of the requirements for the degree of \\
		Doctor of Philosophy \\
		in Physics \\
		2016
      \end{center}
      \thispagestyle{empty}
     \clearpage
     \ncsu@defaultspacing
}

\newcommand{\makecovertwo}{
     \clearpage
     \singlespacing
      \begin{center}
        \begin{Large}
        \textbf{\ncsu@thesistitle} \\[2em]
        \end{Large}
               
        by \\ [2em]
        \@author \\
        M.S., Marquette University, 2011 \\
        B.S., Sichuan University, 2009 \\ [2em]
        Advisor: Stefan Boettcher, Ph.D. \\[5em]
        A dissertation submitted to the Faculty of the \\
        James T. Laney School of Graduate Studies of Emory University \\
        in partial fulfillment of the requirements for the degree of \\
		Doctor of Philosophy \\
		in Physics \\
		2016
      \end{center}
       \thispagestyle{empty}
     \clearpage
     \ncsu@defaultspacing
}

%% -------------------------------------------------------------------------- %%
%% ------------ Dedication, biography, and acknowledgements ----------------- %%
%% -------------------------------------------------------------------------- %%
%%  Template
\newcommand{\ncsu@dedbioack}[1]{%
    \chapter*{#1}%
}

%%  Biography
\newcommand{\biographytitle}{Biography}
\newenvironment{biography}
{\ncsu@dedbioack{\biographytitle}}
{\clearpage}

%%  Acknowledgements
\newcommand{\acknowledgementstitle}{Acknowledgements}
\newenvironment{acknowledgements}
{\ncsu@dedbioack{\acknowledgementstitle}}
{\thispagestyle{empty} \clearpage}

%% Citation to previous work
\newcommand{\citationtitle}{Citations to Previously Published Work}
\newenvironment{prev_citation}{
\singlespacing
{\ncsu@dedbioack{\citationtitle}}
}
{\thispagestyle{empty} \clearpage}

%% -------------------------------------------------------------------------- %%
%% -------------------------  Table of Contents ----------------------------- %%
%% -------------------------------------------------------------------------- %%

\renewcommand{\contentsname}{Table of Contents}
\settowidth{\cftchapnumwidth}{\textbf{Chapter 00}}
\renewcommand{\cftchappresnum}{Chapter }
\renewcommand{\cftchapleader}{\cftdotfill{\cftsecdotsep}}
\newcommand{\thesistableofcontents}{%
    \singlespacing%
    \setlength{\cftbeforetoctitleskip}{-\baselineskip}
    \setlength{\ncsu@leading}{\baselineskip-\ncsu@capitalheight}%
    \setlength{\cftaftertoctitleskip}%
              {\ncsu@afterschapsep-\cftbeforechapskip-\ncsu@leading}%
              
    \addtocontents{toc}{\protect\thispagestyle{empty}}
	\tableofcontents
	\thispagestyle{empty}
  	\clearpage
    \ncsu@defaultspacing%
}
%% tocloft typesets ToC title, so we must give it our heading formats
\renewcommand{\cfttoctitlefont}{\hfill\ncsu@headingformat}
\renewcommand{\cftaftertoctitle}{\hfill}

\renewcommand{\listfigurename}{List of Figures}
%\renewcommand{\listfigurename}{LIST OF FIGURES} %SPB
\settowidth{\cftfignumwidth}{Figure 00.00}

\renewcommand{\cftfigpresnum}{Figure }
\newcommand{\thesislistoffigures}{%
    \singlespacing%
    \setlength{\cftbeforeloftitleskip}{-\baselineskip}
    \setlength{\cftfigindent}{.01in} %SPB
    \setlength{\ncsu@leading}{\baselineskip-\ncsu@capitalheight}%
    \setlength{\cftafterloftitleskip}%
              {\ncsu@afterschapsep-\cftbeforechapskip-\ncsu@leading}%
              
    \cleardoublepage
	\begingroup
	\makeatletter
	\let\ps@plain\ps@empty
	\makeatother
	\pagestyle{empty}
  	\listoffigures %
  	\cleardoublepage
	\endgroup

    \ncsu@defaultspacing%
}
\renewcommand{\cftloftitlefont}{\hfill\ncsu@headingformat}
\renewcommand{\cftafterloftitle}{\hfill}

\renewcommand{\listtablename}{List of Tables} %SPB
%\renewcommand{\listtablename}{LIST OF TABLES} %SPB
\settowidth{\cfttabnumwidth}{Table 00.00}
\renewcommand{\cfttabpresnum}{Table }
\newcommand{\thesislistoftables}{%
    \singlespacing%
    \setlength{\cftbeforelottitleskip}{-\baselineskip}
    \setlength{\cfttabindent}{.01in} %SPB
    \setlength{\ncsu@leading}{\baselineskip-\ncsu@capitalheight}%
    \setlength{\cftafterlottitleskip}%
              {\ncsu@afterschapsep-\cftbeforechapskip-\ncsu@leading}%
              
    \cleardoublepage
	\begingroup
	\makeatletter
	\let\ps@plain\ps@empty
	\makeatother
	\pagestyle{empty}
    \listoftables%
    \cleardoublepage
	\endgroup

    \clearpage%
    \ncsu@defaultspacing%
}
\renewcommand{\cftlottitlefont}{\hfill\ncsu@headingformat}
\renewcommand{\cftafterlottitle}{\hfill}

%% -------------------------------------------------------------------------- %%
%% ------------------------- Landscape Pages  ------------------------------- %%
%% -------------------------------------------------------------------------- %%
%%  Loading pdflscape package for page rotation
\RequirePackage{pdflscape}
%%  Loading fancyhdr package to change page numbering
\RequirePackage{fancyhdr}
%%  Loading sideways for sideways environment
\RequirePackage{rotating}

\newlength{\portraittextheight}
\setlength{\portraittextheight}{\textheight}
\addtolength{\portraittextheight}{\footskip}
\newlength{\portraittextwidth}
\setlength{\portraittextwidth}{\textwidth}

\RequirePackage{afterpage}
\RequirePackage{changepage}

%% Rotate page and page number
\newenvironment{lscapenum}
{
%%   The standard landscape environment does not rotate the margins and
%%   page number.  However, we want the margins (mostly) rotated as well
%%   as the page number.  To change the margins, we use the \changepage
%%   command from the \changepage package.  Since the \textheight and
%%   \textwidth dimensions get switched, we need to reduce the rotated
%%   \textheight by \ncsu@footskip and increase the \textwidth by the
%%   same measure.
%%
%%   The syntax of the \changepage command is
%% \changepage{textheight}{textwidth}{evensidemargin}{oddsidemargin}{columnsep}
%%            {topmargin}{headheight}{headsep}{footskip}
   \clearpage
   \changepage{\ncsu@footskip}{-\ncsu@footskip}{}{}{}{}{}{}{}
   \begin{landscape}
     \pagestyle{fancy}
     \fancyhead{} % clear all header fields
     \fancyfoot{} % clear all footer fields
     \renewcommand{\headrulewidth}{0pt}
     \rfoot{
       %% Since the header and footer are unrotated, we must manually
       %% place the page number.  The picture environment is included in
       %% standard LaTeX.
       \setlength{\unitlength}{\portraittextheight}
       \begin{picture}(0,0)
           \put(0.055,0.557){\begin{sideways}\thepage\end{sideways}} % page number
       \end{picture}
     }

}
{
   \end{landscape}
   %% The \textheight is too large leaving landscape environment
   %%  since we invoked the \changepage.  We manually reduce it.
   \addtolength{\textheight}{-\ncsu@footskip}

   %% \textheight controls 3 other variables, which must be updated.
   %% See \endlandscape definition in lscape.sty
   \global\@colht=\textheight
   \global\vsize=\textheight
   \global\@colroom=\textheight

   \pagestyle{plain}
}

%% Rotate page but NOT page number
\newenvironment{lscapenotnum}
{
   \begin{landscape}
}
{
   \end{landscape}
}

%% -------------------------------------------------------------------------- %%
%% ---------------------------  Captions      ------------------------------- %%
%% -------------------------------------------------------------------------- %%
%%  Adding space to start of caption
%%   save
\RequirePackage{caption}
\DeclareCaptionLabelSeparator{colon}{:\hspace{0.8em}}

%% -------------------------------------------------------------------------- %%
%% -------------------------  References ------------------------------------ %%
%% -------------------------------------------------------------------------- %%
\renewcommand{\bibname}{References}
%% natbib redefines \thebibliography, so we add our own requirements
%% at \begin{document} time.
\AtBeginDocument{
    \let\ncsu@thebibliography\thebibliography
    \let\endncsu@thebibliography\endthebibliography
    \renewenvironment{thebibliography}[1]
    {
        \doublespacing %%modification 12/12/13 by isbayram@ncsu.edu, Grad School requires double spacing in the references
        \begin{ncsu@thebibliography}{#1}
    }
    {
        \end{ncsu@thebibliography}
        \ncsu@defaultspacing
    }
}

%% -------------------------------------------------------------------------- %%
%% ------------------------- The Appendices  -------------------------------- %%
%% -------------------------------------------------------------------------- %%

\newlength{\appendixwidth}
\settowidth{\appendixwidth}{Appendix A}
\addtolength{\appendixwidth}{1em}
\providecommand{\appendixheader}[1]{
    \renewcommand{\ncsu@appendixheader}{#1}
}

%% Count the total number of appendices
\RequirePackage{totcount}
\regtotcounter{chapter} % chapter counter is used for appendices

\AtBeginDocument{%
%% Since ETD requirements are different than LaTeX's normal formatting of
%% the Table of Content's Appendix listings,
%% we change the indentation of chapters, sections, and subsections in
%% the appendix part of the table of contents, but no deeper divisions,
%% so we must cap the tocdepth at 2:
\newcounter{ncsu@appendixtocdepth}
\ifthenelse{\c@tocdepth > 2}
{\setcounter{ncsu@appendixtocdepth}{2}}
{\setcounter{ncsu@appendixtocdepth}{\c@tocdepth}}

\let\ncsu@appendix\appendix
\renewcommand{\appendix}{
    \ifthenelse{\isundefined{\ncsu@appendixheader}}
    % Appendix title page heading not defined, define it
    {
        \ifthenelse{\totvalue{chapter} = 1 }%
        %% number of appendices = 1
        {\newcommand{\ncsu@appendixheader}{Appendix}}
        %% number of appendices > 1
        {\newcommand{\ncsu@appendixheader}{Appendices}}
    }
    % Appendix title page heading already defined, do nothing
    {}

    \chapter*{\ncsu@appendixheader} %% insert 'Appendix' titlepage
    \addcontentsline{toc}{chapter}{\ncsu@appendixheader}
    %% Adding "Appendix" before letter in TOC
    %%  We must change the width of the box which the appendix letter is set in
    %%  and add the "Appendix " before the letter.  See tocloft documentation.
    \addtocontents{toc}{
        \protect\renewcommand{\protect\cftchappresnum}{Appendix }
        \protect\renewcommand{\protect\cftchapnumwidth}{\protect\appendixwidth}
        \protect\renewcommand{\protect\cftchapfont}{\protect\cftsecfont}
        \protect\setlength{\protect\cftbeforechapskip}{\protect\cftbeforesecskip}
        \protect\setlength{\protect\cftchapindent}{\protect\cftsecindent}
        \protect\renewcommand{\protect\cftchappagefont}{\protect\cftchapfont}
        \protect\setlength{\protect\cftsecindent}{\protect\cftsubsecindent}
        \protect\addtolength{\protect\cftsubsecindent}{\cftsecnumwidth}
        \protect\setcounter{tocdepth}{\c@ncsu@appendixtocdepth}
    }
    \clearpage
    \ncsu@appendix %% call internal \appendix macro
} % \renewcommand{appendix}
} % \AtBeginDocument

\ncsu@defaultspacing
%% 
%% Copyright (C) 2009-2012 by Josh Hykes <jmhykes@ncsu.edu> and others
%% 
%% This work may be distributed and/or modified under the
%% conditions of the LaTeX Project Public License (LPPL), either
%% version 1.3c of this license or (at your option) any later
%% version.  The latest version of this license is in the file:
%% 
%% http://www.latex-project.org/lppl.txt
%% 
%% This work is "maintained" (as per LPPL maintenance status).
%% 
%% This work consists of the file  ncsuthesis.dtx
%% and the derived files           ncsuthesis.ins,
%%                                 ncsuthesis.pdf and
%%                                 ncsuthesis.cls.
%% 
%%
%% End of file `ncsuthesis.cls'.
