TARGET_HTML  = ABMPTools-user-manual.html
TARGET_PDF   = ABMPTools-user-manual.pdf
MARKDOWN     = ABMPTools-user-manual.md
FLAGS_HTML   = -c github.css -t html5 -s --self-contained
FLAGS_PDF    = -c github.css -s --pdf-engine=lualatex -H h-lulatexja.tex

all: $(TARGET_HTML) #$(TARGET_PDF)
	#

$(TARGET_HTML):$(MARKDOWN)
	pandoc -f gfm  $(MARKDOWN) $(FLAGS_HTML) -o $(TARGET_HTML)

$(TARGET_PDF):$(MARKDOWN)
	pandoc -f gfm  $(MARKDOWN) $(FLAGS_PDF)  -o $(TARGET_PDF)

clean:
	rm $(TARGET_HTML) $(TARGET_PDF)
