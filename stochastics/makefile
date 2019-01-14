CC=gcc
CFLAGS=-lgsl -lcblas -lm

stochastics: adaptivewaves_changebinaryfile.c  adaptivewaves_texttobinary.c  travelingwavepeak_diff_event.c  travelingwavepeak_exp_event.c adaptivewaves_readfile.c          travelingwavepeak_diff.c      travelingwavepeak_exp.c
	$(CC) $(CFLAGS) adaptivewaves_changebinaryfile.c   -o adaptivewaves_changebinaryfile
	$(CC) $(CFLAGS) adaptivewaves_readfile.c           -o adaptivewaves_readfile
	$(CC) $(CFLAGS) adaptivewaves_texttobinary.c       -o adaptivewaves_texttobinary
	$(CC) $(CFLAGS) travelingwavepeak_diff.c           -o travelingwavepeak_diff
	$(CC) $(CFLAGS) travelingwavepeak_diff_event.c     -o travelingwavepeak_diff_event
	$(CC) $(CFLAGS) travelingwavepeak_exp.c            -o travelingwavepeak_exp
	$(CC) $(CFLAGS) travelingwavepeak_exp_event.c      -o travelingwavepeak_exp_event

clean:
	rm -rf adaptivewaves_changebinaryfile adaptivewaves_readfile adaptivewaves_texttobinary
	rm -rf travelingwavepeak_diff travelingwavepeak_diff_event travelingwavepeak_exp travelingwavepeak_exp_event

