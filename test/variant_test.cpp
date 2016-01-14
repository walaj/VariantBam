#include <climits>
#include <boost/test/unit_test.hpp>

#include "SnowTools/MiniRules.h"
#include "VariantBamWalker.h"

BOOST_AUTO_TEST_CASE( example_case_1 ) {

  SnowTools::BamWalker bw("small.bam");

  std::string rules = "mlregion@test.vcf%all";
  SnowTools::MiniRulesCollection mr(rules, bw.header());

  BOOST_CHECK_EQUAL( mr.numRules(), 1 );
  BOOST_CHECK_EQUAL( mr.m_regions[0].pad, 0);
  BOOST_TEST( !mr.m_regions[0].excluder );
  BOOST_TEST( mr.m_regions[0].m_applies_to_mate );
  BOOST_CHECK_EQUAL( mr.m_regions[0].m_grv.size(), 3);
  BOOST_CHECK_EQUAL( mr.m_regions[0].mrc, &mr);

}

BOOST_AUTO_TEST_CASE( example_case_2 ) {

  SnowTools::BamWalker bw("small.bam");
  std::string rules = "region@WG%phred[4,100];length[50,1000];mapq[1,60];!duplicate;!hardclip;!qcfail";
  SnowTools::MiniRulesCollection mr(rules, bw.header());

  BOOST_CHECK_EQUAL( mr.numRules(), 1 );
  BOOST_CHECK_EQUAL( mr.m_regions[0].pad, 0);
  BOOST_TEST( !mr.m_regions[0].excluder );
  BOOST_TEST( !mr.m_regions[0].m_applies_to_mate );
  BOOST_TEST( mr.m_regions[0].m_whole_genome );
  BOOST_CHECK_EQUAL( mr.m_regions[0].m_grv.size(), 0);

  //BOOST_CHECK_EQUAL( mr.m_regions[0].mrc, &mr);
  BOOST_CHECK_EQUAL( mr.m_regions[0].m_abstract_rules.size(), 1);  
  BOOST_CHECK_EQUAL( mr.m_regions[0].m_abstract_rules[0].phred.min, 4);  
  BOOST_CHECK_EQUAL( mr.m_regions[0].m_abstract_rules[0].phred.max, 100);  
  BOOST_CHECK_EQUAL( mr.m_regions[0].m_abstract_rules[0].len.min, 50);  
  BOOST_CHECK_EQUAL( mr.m_regions[0].m_abstract_rules[0].len.max, 1000);  
  BOOST_CHECK_EQUAL( mr.m_regions[0].m_abstract_rules[0].mapq.min, 1);  
  BOOST_CHECK_EQUAL( mr.m_regions[0].m_abstract_rules[0].mapq.max, 60);  

  // test the flags
  BOOST_TEST( mr.m_regions[0].m_abstract_rules[0].fr.dup.isOff());  
  BOOST_TEST( !mr.m_regions[0].m_abstract_rules[0].fr.dup.isOn());  
  BOOST_TEST( !mr.m_regions[0].m_abstract_rules[0].fr.dup.isNA());  
  BOOST_TEST( mr.m_regions[0].m_abstract_rules[0].fr.hardclip.isOff());  
  BOOST_TEST( !mr.m_regions[0].m_abstract_rules[0].fr.hardclip.isOn());  
  BOOST_TEST( !mr.m_regions[0].m_abstract_rules[0].fr.hardclip.isNA());  
  BOOST_TEST( mr.m_regions[0].m_abstract_rules[0].fr.qcfail.isOff());  
  BOOST_TEST( !mr.m_regions[0].m_abstract_rules[0].fr.qcfail.isOn());  
  BOOST_TEST( !mr.m_regions[0].m_abstract_rules[0].fr.qcfail.isNA());  

} 

BOOST_AUTO_TEST_CASE( example_case_3 ) {

  SnowTools::BamWalker bw("small.bam");  
  std::string rules = "global@nbases[0,0];!hardclip;!supplementary;!duplicate;!qcfail;phred[4,100];%region@WG%discordant[0,1000];mapq[1,1000]%mapq[1,1000];clip[5,1000]%ins[1,1000];mapq[1,100]%del[1,1000];mapq[1,1000]";
  SnowTools::MiniRulesCollection mr(rules, bw.header());

  // discordant rule generates 5 total rules 
  BOOST_CHECK_EQUAL( mr.numRules(), 8 ); //(3 + 5)
  BOOST_CHECK_EQUAL( mr.m_regions.size(), 1 );
  BOOST_CHECK_EQUAL( mr.m_regions[0].size(), 8 ); // size is m_abstract_rules size
  BOOST_CHECK_EQUAL( mr.m_regions[0].pad, 0);
  BOOST_TEST( !mr.m_regions[0].excluder );
  BOOST_TEST( !mr.m_regions[0].m_applies_to_mate );
  BOOST_TEST( mr.m_regions[0].m_whole_genome );
  BOOST_CHECK_EQUAL( mr.m_regions[0].m_grv.size(), 0);

  // check the discordant rules
  SnowTools::AbstractRule * ar = &(mr.m_regions[0].m_abstract_rules[0]);
  BOOST_CHECK_EQUAL( ar->isize.min, 0 ); 
  BOOST_CHECK_EQUAL( ar->isize.max, 1000 ); 
  BOOST_TEST( ar->isize.inverted ); 
  ar = &(mr.m_regions[0].m_abstract_rules[1]);
  BOOST_TEST(  ar->fr.ff.isOn() ); 
  BOOST_TEST( !ar->fr.ff.isOff() ); 
  BOOST_TEST( !ar->fr.ff.isNA() ); 

  BOOST_CHECK_EQUAL( mr.m_regions[0].m_abstract_rules[0].phred.min, 4);  
  BOOST_CHECK_EQUAL( mr.m_regions[0].m_abstract_rules[0].phred.max, 100);  
  BOOST_CHECK_EQUAL( mr.m_regions[0].m_abstract_rules[0].mapq.min, 1);  
  BOOST_CHECK_EQUAL( mr.m_regions[0].m_abstract_rules[0].mapq.max, 1000);  

  // test the flags
  SnowTools::FlagRule * fr = &(mr.m_regions[0].m_abstract_rules[0].fr);
  BOOST_TEST(  fr->dup.isOff());  // dup is OFF
  BOOST_TEST( !fr->dup.isOn());  
  BOOST_TEST( !fr->dup.isNA()); 
  BOOST_TEST(  fr->supp.isOff());  // supp is OFF
  BOOST_TEST( !fr->supp.isOn());  
  BOOST_TEST( !fr->supp.isNA());  
  BOOST_TEST(  fr->mapped.isNA());  // mapped is NA
  BOOST_TEST( !fr->mapped.isOn());  
  BOOST_TEST( !fr->mapped.isOff());  
  BOOST_TEST(  fr->hardclip.isOff());  // hardclip is OFF
  BOOST_TEST( !fr->hardclip.isOn());  
  BOOST_TEST( !fr->hardclip.isNA());  
  BOOST_TEST(  fr->qcfail.isOff()); // qcfail is OFF 
  BOOST_TEST( !fr->qcfail.isOn());  
  BOOST_TEST( !fr->qcfail.isNA());  

}
