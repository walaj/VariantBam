#include "VariantBamWalker.h"

using SnowTools::ReadCount;

void VariantBamWalker::writeVariantBam() 
{

#ifndef __APPLE__
  // start the timer
  clock_gettime(CLOCK_MONOTONIC, &start);
#endif

  ReadCount rc_main;
  SnowTools::BamRead r;
  bool rule;

  while (GetNextRead(r, rule))
    {

      //std::cerr << "...r " << r << " rule " << rule << std::endl;
      TrackSeenRead(r);
      // read is valid
      if (rule) {
	++rc_main.keep;

	if (fop) // if we specified an output file, write it
	  WriteAlignment(r);
      } 

      if (++rc_main.total % 1000000 == 0 && m_verbose)
	printMessage(rc_main, r);

    }
  if (m_verbose)
    printMessage(rc_main,r);
}

void VariantBamWalker::TrackSeenRead(SnowTools::BamRead &r)
{
  m_stats.addRead(r);
}

void VariantBamWalker::printMessage(const ReadCount &rc_main, const SnowTools::BamRead &r) const 
{
  
  char buffer[120];
  std::string posstring = SnowTools::AddCommas<int>(r.Position());
  std::sprintf (buffer, "Reading read %11s at %2s:%-11s. Kept %10s (%2d%%) [running count]",  
		rc_main.totalString().c_str(), SnowTools::GenomicRegion::chrToString(r.ChrID()).c_str(), posstring.c_str(),  
		rc_main.keepString().c_str(), rc_main.percent());
  //buffer[99] = '\0';
  std::string outr(buffer);
  //std::printf ("%s | ",buffer);
  std::cerr << outr << SnowTools::displayRuntime(start) << std::endl;;
  
}
