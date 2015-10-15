#include <iostream>
#include <map>
#include <string>

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"

int FindMin( int n, Float_t *a);
int FindMax( int n, Float_t *a);
float LinearFit_Baseline(TH1F * pulse, const int index_min, const int range);


int main(int argc, char *argv[]){

  //get the input file
  TFile *file = file = new TFile(argv[1]);  
  TTree* H4tree = (TTree*)file->Get("H4tree");

  // get the tree variables
  int ndigisamples = 18432;
  Float_t cv[18432];
  UInt_t ci[18432];
  UInt_t cgroup[18432];
  UInt_t cchannel[18432];

  H4tree->SetBranchAddress("digiSampleValue", cv);
  H4tree->SetBranchAddress("digiSampleIndex", ci);
  H4tree->SetBranchAddress("digiGroup", cgroup);
  H4tree->SetBranchAddress("digiChannel", cchannel);

  //define some histograms
  float LinearFit_Baseline(TH1F * pulse, const int index_min, const int range);
  TH1F *CH1pulse   = new TH1F("CH1pulse","CH1pulse",1024,0,1024);
  TH1F *CH2pulse   = new TH1F("CH2pulse","CH2pulse",1024,0,1024);
  TH1F *CH3pulse   = new TH1F("CH3pulse","CH3pulse",1024,0,1024);
  TH1F *CH4pulse   = new TH1F("CH4pulse","CH4pulse",1024,0,1024);
  TH1F *CH5pulse   = new TH1F("CH5pulse","CH5pulse",1024,0,1024);
  TH1F *CH6pulse   = new TH1F("CH6pulse","CH6pulse",1024,0,1024);
  TH1F *CH7pulse   = new TH1F("CH7pulse","CH7pulse",1024,0,1024);
  TH1F *CH8pulse   = new TH1F("CH8pulse","CH8pulse",1024,0,1024);

  // Create the output file with a TTree
  float ch1Amp = 0;
  float ch2Amp = 0;
  float ch3Amp = 0;
  float ch4Amp = 0;
  float ch5Amp = 0;
  float ch6Amp = 0;
  float ch7Amp = 0;
  float ch8Amp = 0;

  TFile fout(argv[2],"recreate");
  TTree *treeOut = new TTree("tree","tree");
  treeOut->Branch("ch1Amp",&ch1Amp,"ch1Amp/F");
  treeOut->Branch("ch2Amp",&ch2Amp,"ch2Amp/F");
  treeOut->Branch("ch3Amp",&ch3Amp,"ch3Amp/F");
  treeOut->Branch("ch4Amp",&ch4Amp,"ch4Amp/F");
  treeOut->Branch("ch5Amp",&ch5Amp,"ch5Amp/F");
  treeOut->Branch("ch6Amp",&ch6Amp,"ch6Amp/F");
  treeOut->Branch("ch7Amp",&ch7Amp,"ch7Amp/F");
  treeOut->Branch("ch8Amp",&ch8Amp,"ch8Amp/F");

  // Get NEvents
  Long64_t nevents = H4tree->GetEntries();

  std::cout<<"Number of events: "<<nevents<<std::endl;

  for(Long64_t i = 0; i<nevents; i++)
    {

      if(i%100==0) std::cout<<"Processing Event: "<<i<<" out of: "<<nevents<<std::endl;
      H4tree->GetEntry(i);
      
      // create the arrays to store ADC counts per channel
      Float_t Channel1Voltages_[1024];
      Float_t Channel2Voltages_[1024];
      Float_t Channel3Voltages_[1024];
      Float_t Channel4Voltages_[1024];
      Float_t Channel5Voltages_[1024];
      Float_t Channel6Voltages_[1024];
      Float_t Channel7Voltages_[1024];
      Float_t Channel8Voltages_[1024];

      // initiallize them
      for(int kk=0;kk<1024;kk++)
	{
	  Channel1Voltages_[kk] = 0.;
	  Channel2Voltages_[kk] = 0.;
	  Channel3Voltages_[kk] = 0.;
	  Channel4Voltages_[kk] = 0.;
	  Channel5Voltages_[kk] = 0.;
	  Channel6Voltages_[kk] = 0.;
	  Channel7Voltages_[kk] = 0.;
	  Channel8Voltages_[kk] = 0.;
	}

      // get the pulses (AKA as digis) from the input tree
      for(int jj =0; jj<ndigisamples; jj++)
	{
	  if(cgroup[jj]==0)
	    {
	      if(cchannel[jj]==1)
		Channel1Voltages_[ci[jj]] = cv[jj];
	      
	      if(cchannel[jj]==2)
		Channel2Voltages_[ci[jj]] = cv[jj];

	      if(cchannel[jj]==3)
		Channel3Voltages_[ci[jj]] = cv[jj];

	      if(cchannel[jj]==4)
		Channel4Voltages_[ci[jj]] = cv[jj];

	      if(cchannel[jj]==5)
		Channel5Voltages_[ci[jj]] = cv[jj];

	      if(cchannel[jj]==6)
		Channel6Voltages_[ci[jj]] = cv[jj];

	      if(cchannel[jj]==7)
		Channel7Voltages_[ci[jj]] = cv[jj];

	      if(cchannel[jj]==8)
		Channel8Voltages_[ci[jj]] = cv[jj];
	    }
	}

      // Fill the pulse histograms, to get the baseline
      for (int ii=0;ii<1024;ii++)
	{
	  CH1pulse->SetBinContent(ii+1,Channel1Voltages_[ii]);
	  CH2pulse->SetBinContent(ii+1,Channel2Voltages_[ii]);
	  CH3pulse->SetBinContent(ii+1,Channel3Voltages_[ii]);
	  CH4pulse->SetBinContent(ii+1,Channel4Voltages_[ii]);
	  CH5pulse->SetBinContent(ii+1,Channel5Voltages_[ii]);
	  CH6pulse->SetBinContent(ii+1,Channel6Voltages_[ii]);
	  CH7pulse->SetBinContent(ii+1,Channel7Voltages_[ii]);
	  CH8pulse->SetBinContent(ii+1,Channel8Voltages_[ii]);
	}

      // get the baseline, fit is done on bins 5 to 10
      float base1 = LinearFit_Baseline( CH1pulse, 5, 10 );
      float base2 = LinearFit_Baseline( CH2pulse, 5, 10 );
      float base3 = LinearFit_Baseline( CH3pulse, 5, 10 );
      float base4 = LinearFit_Baseline( CH4pulse, 5, 10 );
      float base5 = LinearFit_Baseline( CH5pulse, 5, 10 );
      float base6 = LinearFit_Baseline( CH6pulse, 5, 10 );
      float base7 = LinearFit_Baseline( CH7pulse, 5, 10 );
      float base8 = LinearFit_Baseline( CH8pulse, 5, 10 );
      
      // subtract the baseline
      for(int ll=0;ll<1024;ll++)
	{
	  Channel1Voltages_[ll] = Channel1Voltages_[ll] - base1;
	  Channel2Voltages_[ll] = Channel2Voltages_[ll] - base2;
	  Channel3Voltages_[ll] = Channel3Voltages_[ll] - base3;
	  Channel4Voltages_[ll] = Channel4Voltages_[ll] - base4;
	  Channel5Voltages_[ll] = Channel5Voltages_[ll] - base5;
	  Channel6Voltages_[ll] = Channel6Voltages_[ll] - base6;
	  Channel7Voltages_[ll] = Channel7Voltages_[ll] - base7;
	  Channel8Voltages_[ll] = Channel8Voltages_[ll] - base8;
	}

      // find the bin with minimum voltage
      int index_min1 = FindMin(1024, Channel1Voltages_); // return index of the min
      int index_min2 = FindMin(1024, Channel2Voltages_);  
      int index_min3 = FindMin(1024, Channel3Voltages_);  
      int index_min4 = FindMin(1024, Channel4Voltages_);  
      int index_min5 = FindMin(1024, Channel5Voltages_);  
      int index_min6 = FindMin(1024, Channel6Voltages_);  
      int index_min7 = FindMin(1024, Channel7Voltages_);  
      int index_min8 = FindMin(1024, Channel8Voltages_);  
      
      // assign the values to tree variables
      ch1Amp = -1 * Channel1Voltages_[index_min1];
      ch2Amp = -1 * Channel2Voltages_[index_min2];
      ch3Amp = -1 * Channel3Voltages_[index_min3];
      ch4Amp = -1 * Channel4Voltages_[index_min4];
      ch5Amp = -1 * Channel5Voltages_[index_min5];
      ch6Amp = -1 * Channel6Voltages_[index_min6];
      ch7Amp = -1 * Channel7Voltages_[index_min7];
      ch8Amp = -1 * Channel8Voltages_[index_min8];

      //Fill the tree
      treeOut->Fill();      
    }

  treeOut->Write();
  file->Close();
}

int FindMin( int n, Float_t *a){
  if (n <= 0 || !a) return -1;
  float xmin = a[5];
  int loc = 0;
  for  (int i = 5; i < n-5; i++) {
    if (xmin > a[i] && a[i+1] < 0.5*a[i])  {
      xmin = a[i];
      loc = i;

    }
  }
  
  return loc;
}

float LinearFit_Baseline(TH1F * pulse, const int bin1, const int bin2)
{
  TF1 *fBaseline = new TF1("fBaseline","pol0",bin1, bin2);
  pulse->Fit("fBaseline","Q","", bin1, bin2);
  float base = fBaseline->GetParameter(0);
  
  return base;
}
