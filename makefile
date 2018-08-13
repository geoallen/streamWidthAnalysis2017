all:
	./smallStreamsAnalysis.r

# fig1_widthMap.pdf: fig1_widthMap.R
#

clean:
	rm -rf figures/
	rm -rf tables/

.PHONY: all clean
