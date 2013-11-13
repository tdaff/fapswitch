"""

Synthetic accessibility scoring functions based on:

Estimation of Synthetic Accessibility Score of Drug-like Molecules based on
Molecular Complexity and Fragment Contributions
Peter Ertl and Ansgar Schuffenhauer
Journal of Cheminformatics 1:8 (2009)
doi:10.1186/1758-2946-1-8

several small modifications to the original paper are included
particularly slightly different formula for marocyclic penalty
and taking into account also molecule symmetry (fingerprint density)

This implementation is derived from the RDKit script but has been
adapted for fapswitch.

"""


import math
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from fapswitch.config import debug
from fapswitch.extensions.fragments import moldb


def ring_analysis(molecule):
    """
    Return the numbers of bridgeheads, spiro and macrocycles.

    """

    ring_info = molecule.GetRingInfo()

    atom_rings = [set(x) for x in ring_info.AtomRings()]
    spiros = set()

    for idx, atom_ring in enumerate(atom_rings):
        for jdx in range(idx + 1, len(atom_rings)):
            shared = atom_ring & atom_rings[jdx]
            if len(shared) == 1:
                spiros.update(shared)

    num_spiro = len(spiros)

    # find bonds that are shared between rings that share at least 2 bonds:
    bond_rings = [set(x) for x in ring_info.BondRings()]
    bridges = set()
    for idx, bond_ring in enumerate(bond_rings):
        for jdx in range(idx + 1, len(bond_rings)):
            shared = bond_ring & bond_rings[jdx]
            if len(shared) > 1:
                atom_counts = defaultdict(int)
                for bond_id in shared:
                    bond = molecule.GetBondWithIdx(bond_id)
                    atom_counts[bond.GetBeginAtomIdx()] += 1
                    atom_counts[bond.GetEndAtomIdx()] += 1
                for atom_id, count in atom_counts.items():
                    if count == 1:
                        bridges.add(atom_id)

    # Count how many macrocycles there are
    num_macrocycles = 0
    for ring in ring_info.AtomRings():
        if len(ring) > 8:
            num_macrocycles += 1

    return len(bridges), num_spiro, num_macrocycles


def sa_score(smiles):
    """
    Return the SA Score for the given smiles representation.

    """

    molecule = Chem.MolFromSmiles(smiles)

    #
    # fragment score
    #

    # use a radius of 2 for circular fingerprint
    try:
        fingerprint = rdMolDescriptors.GetMorganFingerprint(molecule, 2)
        fingerprint = fingerprint.GetNonzeroElements()
    except Exception as error:
        # Will throw a boost error for N+ so we just give a 0 for score
        debug(error)
        return 0

    fragment_score = 0.0
    fragment_count = 0

    # Count frequencies of fragments
    for bit_id, count in fingerprint.items():
        fragment_count += count
        fragment_score += moldb.get(bit_id, -4) * count

    fragment_score /= fragment_count

    #
    # features score
    #

    num_atoms = molecule.GetNumAtoms()
    num_chiral_centers = len(Chem.FindMolChiralCenters(molecule,
                                                       includeUnassigned=True))
    num_bridgeheads, num_spiro, num_macrocycles = ring_analysis(molecule)

    size_penalty = (num_atoms ** 1.005) - num_atoms
    stereo_penalty = math.log10(num_chiral_centers + 1)
    spiro_penalty = math.log10(num_spiro + 1)
    bridge_penalty = math.log10(num_bridgeheads + 1)

    macrocycle_penalty = 0.0
    # ---------------------------------------
    # This differs from the paper, which defines:
    #  macrocycle_penalty = math.log10(num_macrocycles + 1)
    # This form generates better results when 2 or more macrocycles are present
    if num_macrocycles > 0:
        macrocycle_penalty = math.log10(2)

    feature_penalty = (0.0 - size_penalty - stereo_penalty - spiro_penalty -
                       bridge_penalty - macrocycle_penalty)

    #
    # Correction for the fingerprint density.
    # Not in the original publication, added in version 1.1
    # to make highly symmetrical molecules easier to synthetise.
    #
    if num_atoms > len(fingerprint):
        fingerprint_density = math.log(float(num_atoms)/len(fingerprint)) * 0.5
    else:
        fingerprint_density = 0.0

    #
    # Total score
    #
    total_score = fragment_score + feature_penalty + fingerprint_density

    # Transform "raw" value into scale between 1 and 10.
    sa_min = -4.0
    sa_max = 2.5
    total_score = 11.0 - (total_score - sa_min + 1) / (sa_max - sa_min) * 9.0
    # smooth the 10-end
    if total_score > 8.0:
        total_score = 8.0 + math.log(total_score + 1.0 - 9.0)

    if total_score > 10.0:
        total_score = 10.0
    elif total_score < 1.0:
        total_score = 1.0

    return total_score


#
#  Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
