#include "VariantBamWalker.h"

using SnowTools::ReadCount;

void VariantBamWalker::writeVariantBam() 
{

#ifndef __APPLE__
  // start the timer
  clock_gettime(CLOCK_MONOTONIC, &start);
#endif

  ReadCount rc_main;
  Read r;
  bool rule;

  while (GetNextRead(r, rule))
    {
      
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
}

void VariantBamWalker::TrackSeenRead(Read &r)
{
  m_stats.addRead(r);
}

void VariantBamWalker::printMessage(const ReadCount &rc_main, const Read &r) const 
{
  
  char buffer[100];
  std::string posstring = SnowTools::AddCommas<int>(r_pos(r));
  std::sprintf (buffer, "Reading read %11s at position %2s:%-11s. Kept %11s (%2d%%) [running count across whole BAM]",  
		rc_main.totalString().c_str(), SnowTools::GenomicRegion::chrToString(r_id(r)).c_str(), posstring.c_str(),  
		rc_main.keepString().c_str(), rc_main.percent());
  buffer[99] = '\0';
  std::string outr(buffer);
  //std::printf ("%s | ",buffer);
  std::cerr << outr << SnowTools::displayRuntime(start) << std::endl;;
  
}
