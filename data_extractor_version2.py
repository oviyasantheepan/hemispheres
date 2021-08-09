#!/usr/bin/python
from ROOT import TFile
import math


class Jet():
    '''
    This is a "container class" to combine the 4 numbers required to describe a jet into a single "object".
    That makes things a bit nicer to use. 
    '''
    def __init__(self, pT, eta, phi, energy):
        self.pT = pT
        self.eta = eta
        self.phi = phi 
        self.energy = energy

    def __repr__(self):
        return "Jet({0}, {1}, {2}, {3})".format(self.pT, self.eta, self.phi, self.energy)


def main(debug, input_path, event_limit):
    

    # Open the data file 
    input_file = TFile.Open(input_path, "READ")
    tree_name = "trees_SRRPV_"
    tree = input_file.Get(tree_name)

    # Find out how many entries there are in the file
    nentries = tree.GetEntries()

    # Loop over the events in the file. Each iteration of the loop is a new event 
    for ientry in range(nentries):
        tree.GetEntry(ientry)

        print("This is event {0}".format(ientry))

        # Limit number of iterations 
        if (event_limit > 0) and ientry > event_limit: 
            print("Reached maximium number of events. Will exit")
            break

        njets = len(tree.jet_pt)
        print("\tThe event contains {0} jets".format(njets))

        jets = []
        for pT, eta, phi, energy in zip(tree.jet_pt, tree.jet_eta, tree.jet_phi, tree.jet_e):
            a_jet = Jet(pT, eta, phi, energy) 
            jets.append(a_jet)

        thrust_phi = calculate_phi(jets)
        print("\tThe phi of the thrust axis is: {0}".format(thrust_phi))

        # TODO: split the jets into each side, depending on their phi coordinate
        jets_side_A = []
        jets_side_B = []
        for jet in jets:
            if -math.pi*0.5 < jet.phi - thrust_phi < math.pi*0.5:
                jets_side_A.append(jet)
            else:
                jets_side_B.append(jet)
            phi_new = rotate(jets, thrust_phi)

def calculate_phi(jets):

    # TODO: Insert your code to calculate phi based on the jets  
    phi = math.pi
    numerator = 0
    denominator = 0
    for jet in jets:
        px = (jet.pT)*math.cos(jet.phi)
        py = (jet.pT)*math.sin(jet.phi)
        numerator = numerator - 2*px*py
        denominator = denominator + py**2 - px**2
    phi = 0.5*(math.atan(numerator/denominator))
    if math.cos(2*phi)*denominator < -math.sin(2*phi)*numerator:
        k = 1
    else:
        k = 0
    phi = phi + k*math.pi*0.5
    return phi

def rotate(jets, thrust_phi):
    # Rotate hemispheres so normal is horizontal
    phi_new = jet.phi - thrust.phi

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--debug", help="Turn on debug messages", action="store_true", default=False)
    parser.add_argument("-i", "--input_path", help="Path to input file", default="")
    parser.add_argument("-e", "--event_limit", help="Maximum number of events to consider", default=-1)
    args = parser.parse_args()
    debug = args.debug
    input_path = args.input_path
    event_limit = int(args.event_limit)

    if not input_path:
        raise RuntimeError("You must specify the path to the input file")

    main(debug, input_path, event_limit)

