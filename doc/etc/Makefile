
photons.ps:	photons.tex angfun.ps ic_table_CBR.ps
	latex photons
	latex photons	
	dvips -o $@ photons

clean:
	-$(RM) -f .##* *.aux *.log *.dvi
