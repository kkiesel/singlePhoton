LHAINC=-I/user/kiesel/pdfCalculations/rata_scratch/local/include
LHALIB=-L/user/kiesel/pdfCalculations/rata_scratch/local/lib

all:
	g++ -c $(shell root-config --cflags) $(LHAINC) PdfSystematic.cc -o PdfSystematic.o
	g++ $(shell root-config --libs) $(LHALIB) -lLHAPDF PdfSystematic.o -o PdfSystematic

