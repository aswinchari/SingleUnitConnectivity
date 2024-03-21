# Generating networks from post spike and coupling filters from single unit spike trains

The code and data in this repository relate to generating networks from single unit spike trains in multiple steps:
- Single unit spike trains are used to generate post spike filters (PSFs) and coupling filters (CFs) that model the probability of one neuron firing in the next 700ms given that itself (PSF) or another neuron (CF) has just fired. This is done using generalised linear modelling. 
- The normalised dot product of PSF and CF waveforms are used to generate connectivity matrices, which we term _autonomous sychrnogenicity_ and _network synchrogenicity_ respectively.
- The properties of these networks may then be studied to inform us about the role of these neurons in driving the networks towards synchronous activity, as is seen in epilepsy. 

The methods are described in further detail in the manuscript 'Single unit-derived connectivity networks in tuberous sclerosis complex reveal propensity for network hypersynchrony driven by tuber-tuber interactions' by Chari et al. 
