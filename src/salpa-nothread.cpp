// posthocsalpa.cpp

#include "LocalFit.h"
#include "NoiseLevels.h"
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>

constexpr int LOG2BUFSIZE = 16; // 65536 scans

void usage() {
  std::cerr
    << "Usage: posthocsalpa -F samplerate_kHz -c channelcount -C fullcount\n"
    << "                 -t threshold_digi -x threshold_std\n"
    << "                 -l halflength_ms\n"
    << "                 -a asymtime_ms -b blanktime_ms -A Ahead_ms\n"
    << "                 -r rail1_digi[,rail2_digi]\n"
    << "                 -p period_ms -d delay_ms -f forcepeg_ms\n"
    << "                 -P forcepeg_filename\n"
    << "                 -B\n"
    << "\n"
    << "Performs post-hoc artifact filtering using LocalFit.\n"
    << "-F must be given before any of the parameters that are specified in ms\n"
    << "-c specifies the number of channels to process; -C the total number of\n"
    << "   channels in a scan. (For instance, for a UCLA probe, -c might be 128\n"
    << "   whereas -C is 142.)\n"
    << "-t and -x are mutually exclusive.\n"
    << "-n can be used to read a previously recorded noise estimate from disk iff\n"
    << "   -x is used. If -x is given without -n, artifilt will estimate the noise\n"
    << "   anew, based on the first 2 s of the recording, using NoiseLevels.\n"
    << "-l specifies the half-width of the fit window (tau).\n"
    << "-a specifies the size of the beginning of the fit window used for initial\n"
    << "   goodness-of-fit estimation.\n"
    << "-b specifies how much of the first fit is blanked.\n"
    << "-A specifies how many samples ahead a peg is detected.\n"
    << "-r specifies one or two digital values that are treated as rails.\n"
    << "-p, -d, -f can be used to force a peg response if the electrode recording\n"
    << "   does not really hit the rail. Pegs of length f are enforced at t=k*p+d.\n"
    << "-P specifies the name of a file containing sample numbers where pegs should\n"
    << "   be forced. -P and -p/-d are mutually exclusive.\n"
    << "-B enables subtracting of baseline before processing. This is useful\n"
    << "   for numerical stability if baseline is far from zero.\n"
    << "\n"
    << "Default values are:\n"
    << "   F = 30,000, c = C = 64, l = 3 ms,\n"
    << "   a = 0.2 ms, b = 0.4 ms, A = 0.2 ms, x = 3,\n"
    << "   r = -32767,32767, and no forced peg response.\n";
exit(1);
}



class Params {
public:
  int nchans;
  int totalchans;
  int freq_hz;
  int thresh_digi;
  float thresh_std;
  int tau_sams;
  int asym_sams;
  int blank_sams;
  int ahead_sams;
  raw_t rail1, rail2;
  int period_sams;
  int delay_sams;
  int forcepeg_sams;
  char const *forcepeg_filename;
  bool basesub;
public:
  Params() {
    nchans = 0;
    totalchans = 0;
    freq_hz = 25000;
    thresh_digi = 0; // i.e., do not use
    thresh_std = 3;
    tau_sams = 3 * freq_hz / 1000; // 3 ms
    asym_sams = 10;
    blank_sams = 20;
    ahead_sams = 5;
    rail1 = -32767;
    rail2 = 32767;
    period_sams = 0; // i.e., do not use
    delay_sams = 0;
    forcepeg_sams = 0;
    forcepeg_filename = 0;
    basesub = false;
  }
  bool fromArgs(int argc, char **argv) {
    // return true if OK
    while (argc>1) {
      argc--;
      argv++;
      if (argv[0][0]=='-') {
        char letter = argv[0][1];
        char *arg;
        if (argv[0][2]>=32 || letter=='B') {
          arg = argv[0] + 2;
        } else {
          argc--;
          argv++;
          if (argc>0) {
            arg = argv[0];
          } else {
            std::cerr << "Unexpected end of arg list\n";
            return false;
          }
        }
        switch (letter) {
        case 'c': nchans = atoi(arg); break;
        case 'C': totalchans = atoi(arg); break;
        case 't': thresh_digi = atoi(arg); thresh_std = 0; break;
        case 'x': thresh_std = atof(arg); thresh_digi = 0; break;
        case 'F': freq_hz = int(1000*atof(arg)); break;
        case 'l': tau_sams = int(freq_hz * atof(arg) / 1000); break;
        case 'a': asym_sams = int(freq_hz * atof(arg) / 1000); break;
        case 'b': blank_sams = int(freq_hz * atof(arg) / 1000); break;
        case 'A': ahead_sams = int(freq_hz * atof(arg) / 1000); break;
        case 'r': {
          rail1 = atoi(arg);
          char *x = std::strchr(arg, ',');
          rail2 = x ? atoi(x+1) : rail1;
        } break;
        case 'p': period_sams = int(freq_hz * atof(arg) / 1000); break;
        case 'd': delay_sams = int(freq_hz * atof(arg) / 1000); break;
        case 'f': forcepeg_sams = int(freq_hz * atof(arg) / 1000); break;
        case 'P': forcepeg_filename = arg; break;
        case 'B': basesub = true; break;
        default:
          std::cerr << "Unknown parameter: " << letter << "\n";
          return false;
        }
      } else {
        std::cerr << "Unexpected argument: " << *argv << "\n";
        return false;
      }
    }
    if (nchans==0)
      nchans = totalchans;
    else if (totalchans==0)
      totalchans = nchans;
    if (nchans==0)
      nchans = totalchans = 64;
    if (nchans>totalchans)
      return false;
    return true;
  }
};

void crash(char const *x) {
  std::cerr << x << "\n";
  std::exit(2);
}

int main(int argc, char **argv) {
  if ((INFTY + 1) != 0) {
    crash("BUG: Infinity isn't.");
    return 2;
  }
  Params p;
  if (!p.fromArgs(argc, argv)) {
    usage();
    return 1;
  }

  if (p.thresh_digi!=0 && p.thresh_std!=0)
    usage();

  FILE *in = stdin;
  FILE *out = stdout;
  FILE *events = 0;

  if (p.forcepeg_filename) {
    events = std::fopen(p.forcepeg_filename, "r");
    if (!events)
      crash("Cannot open timestamp filename");
  }

  constexpr int BUFSAMS = 1<<LOG2BUFSIZE;
  constexpr int FRAGSAMS = BUFSAMS / 4;
  constexpr int FRAGMASK = FRAGSAMS - 1;
  
  std::vector<raw_t> inbuf(p.totalchans*BUFSAMS);
  std::vector<raw_t> outbuf(p.totalchans*BUFSAMS);  
  std::vector<CyclBuf<raw_t>> inbufs;
  std::vector<CyclBuf<raw_t>> outbufs;
  for (int c=0; c<p.totalchans; c++) {
    inbufs.push_back(CyclBuf<raw_t>(inbuf.data() + c,
                                    LOG2BUFSIZE, p.totalchans));
    outbufs.push_back(CyclBuf<raw_t>(outbuf.data() + c,
                                     LOG2BUFSIZE, p.totalchans));
  }

  timeref_t filledto = 0;
  timeref_t basesubto = 0;
  timeref_t processedto = 0;
  timeref_t savedto = 0;
  timeref_t nextpeg = p.delay_sams ? p.delay_sams : INFTY;

  if (events) 
    if (std::fscanf(events, "%lu", &nextpeg) < 1)
      nextpeg = INFTY;
  
  std::vector<float> thresh(p.nchans, p.thresh_digi);
  std::vector<raw_t> basesub(p.nchans, 0);

  if (p.thresh_std!=0 || p.basesub) {
    int n = std::fread(inbuf.data(),
                       p.totalchans*sizeof(raw_t), FRAGSAMS,
                       in);
    if (n != FRAGSAMS) 
      crash("Cannot read enough data for noise estimate");
    filledto = n;
    for (int c=0; c<p.nchans; c++) {
      NoiseLevels noise;
      noise.train(inbufs[c], 0, n);
      noise.makeready();
      if (p.thresh_std!=0)
        thresh[c] = p.thresh_std * noise.std();
      if (p.basesub)
        basesub[c] = -noise.mean();
    }
  }
  
  std::vector<LocalFit *> fitters(p.nchans, 0);
  for (int c=0; c<p.nchans; c++) {
    fitters[c] = new LocalFit(inbufs[c], outbufs[c],
                              0, thresh[c], p.tau_sams,
                              p.blank_sams, p.ahead_sams,
                              p.asym_sams);
    fitters[c]->setrail(p.rail1 + basesub[c], p.rail2 + basesub[c]);
  }
  
  bool at_eof = false;
  bool go_on = true;

  while (go_on) {
    go_on = false;

    // -- save some stuff
    timeref_t mightsaveto = processedto & ~FRAGMASK;
    while (savedto < mightsaveto) {
      go_on = true;
      fwrite(&(outbufs[0][savedto]),
             sizeof(raw_t)*p.totalchans, FRAGSAMS,
             out);
      savedto += FRAGSAMS;
    }

    // -- subtract baseline
    if (p.basesub) {    
      while (basesubto < filledto) {
        for (int c=0; c<p.nchans; c++)
          inbufs[c][basesubto] += basesub[c];
        basesubto++;
      }
    } else {
      basesubto = filledto;
    }

    // -- subtract artifacts
    timeref_t mightprocessto = filledto
      - p.forcepeg_sams - 3*p.tau_sams - 2;
    if (mightprocessto > savedto + BUFSAMS)
      mightprocessto = savedto + BUFSAMS;
    while (processedto < mightprocessto) {
      go_on=true;
      if (nextpeg < mightprocessto) {
        // process to upcoming peg
        // first, copy non-electrode channels
        for (timeref_t tt=processedto; tt<nextpeg+p.forcepeg_sams; tt++)
          for (int hw=p.nchans; hw<p.totalchans; hw++)
            outbufs[hw][tt] = inbufs[hw][tt];
        // then, process electrode channels
        for (int hw=0; hw<p.nchans; hw++) 
          if (fitters[hw]->forcepeg(nextpeg, nextpeg+p.forcepeg_sams)
              != nextpeg+p.forcepeg_sams)
            crash("LocalFit doesn't like my data (while pegging)!");
        processedto = nextpeg + p.forcepeg_sams;
        // find next peg
        if (p.period_sams) 
          nextpeg += p.period_sams;
        else if (events)
          if (std::fscanf(events, "%lu", &nextpeg) < 1)
            nextpeg = INFTY;
      } else {
        // process as far as we have loaded
        // first, copy non-electrode channels
        for (timeref_t tt=processedto; tt<mightprocessto; tt++)
          for (int hw=p.nchans; hw<p.totalchans; hw++)
            outbufs[hw][tt] = inbufs[hw][tt];
        // then, process electrode channels
        for (int hw=0; hw<p.nchans; hw++)
          if (fitters[hw]->process(mightprocessto) != mightprocessto)
            crash("LocalFit doesn't like my data!");
        processedto=mightprocessto;
      }
    }

    // -- load some data
    if (!at_eof) {
      int n = std::fread(&inbufs[0][filledto],
                         sizeof(raw_t)*p.totalchans, FRAGSAMS,
                         in);
      if (n<0) 
        crash("Cannot read from input");
      filledto += n;
      if (n>0)
        go_on=true;
      if (n != FRAGSAMS)
        at_eof=true;
    }
  }

  // -- EOF!
    
  std::cerr << "go_on\n" << savedto << " " << processedto << " "
            << filledto << " " << nextpeg << " " << events << "\n";
  fitters[0]->report();
  
    // let's process the last bit...
  timeref_t mightprocessto = filledto - 2*p.tau_sams - 1;
  for (timeref_t tt=processedto; tt<mightprocessto; tt++)
    for (int hw=p.nchans; hw<p.totalchans; hw++)
      outbufs[hw][tt] = inbufs[hw][tt];
  for (int hw=0; hw<p.nchans; hw++) 
    if (fitters[hw]->process(mightprocessto) != mightprocessto)
      crash("LocalFit doesn't like my data!");

  // let's save last bit
  for (timeref_t t=savedto; t<processedto; t++)
    std::fwrite(&outbufs[0][t], sizeof(raw_t)*p.totalchans, 1, out);

  return 0;
}

  
      
    
