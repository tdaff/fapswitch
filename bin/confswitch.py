#!/usr/bin/env python

"""
confswitch.py

For a set functionalisation, find which conformer has the lowest energy. Use
faps to calculate the energy.

"""


import json
import math
import os
import random
import re
import string
import subprocess
import sys
from os import path
from os.path import dirname, realpath
sys.path.insert(1, dirname(dirname(realpath(__file__))))

import fapswitch
from fapswitch.backend.cif_file import CifFileBackend
from fapswitch.config import options
from fapswitch.config import debug, info
from fapswitch.core.io import load_structure
from fapswitch.core.methods import site_replace
from fapswitch.functional_groups import functional_groups


def optimise_conformation(structure, site_list, backends):
    """
    Optimise the conformation using a genetic algorithm

    """

    # Set up the files and their names!
    known_name = path.join(os.getcwd(), '{}-known.json'.format(structure.name))
    output_name = path.join(os.getcwd(), '{}-opt.csv'.format(structure.name))
    opt_path = path.join(os.getcwd(), 'opt_{}'.format(structure.name))

    output_file = open(output_name, 'w')
    output_file.write('#generation,most_fit,average\n')

    # Use a library if one exists
    if os.path.exists(known_name):
        known_individuals = json.load(open(known_name))
        info("Existing library found")
    else:
        known_individuals = {}

    # Run in a subdirectory
    try:
        os.mkdir(opt_path)
    except OSError:
        pass
    os.chdir(opt_path)

    elitism = 3
    max_pop = 30
    mate_pop = 15
    mutation_rate = 0.1
    max_generations = 30

    population = []

    generation = 1
    fill_with_randoms(population, max_pop, structure, site_list, backends)

    while True:
        for individual in population:
            info("Fitness of {}".format(individual['gene']))
            if individual['gene'] in known_individuals:
                known_individual = known_individuals[individual['gene']]
                individual['fitness'] = known_individual['fitness']
            else:
                individual['fitness'] = measure_fitness(structure, site_list,
                                                        individual)
            info("From faps {}".format(individual['fitness']))
            known_individuals[individual['gene']] = individual

        population.sort()
        most_fit = population[0]
        average_fit = sum(x['fitness'] for x in population)/len(population)
        info("Most fit {}: {}".format(most_fit['gene'], most_fit['fitness']))
        info("Generation average: {}".format(average_fit))
        output_file.write("{},{},{}\n".format(generation, most_fit['fitness'],
                                              average_fit))
        output_file.flush()

        if generation > max_generations:
            break
        else:
            new_population = population[0:elitism]
            mate_genes(population, new_population, mate_pop, structure,
                       site_list, backends, mutation_rate)
            fill_with_randoms(new_population, max_pop, structure, site_list,
                              backends)
            population = new_population
            generation += 1
            if not generation % 5:
                json.dump(known_individuals, open(known_name, 'w'))

    os.chdir('..')


def mate_genes(old_population, new_population, count, structure, site_list,
               backends, mutation_rate):
    """
    Make some new structures by mating current ones, based on position in
    the population.
    """
    while count > 0:
        index_0 = int(math.floor(len(old_population)*(random.random()**2)))
        index_1 = int(math.floor(len(old_population)*(random.random()**2)))
        old_gene_0 = old_population[index_0]['gene']
        old_gene_1 = old_population[index_1]['gene']
        debug("Mating {} and {}".format(old_gene_0, old_gene_1))
        new_gene = []
        for gene_part_0, gene_part_1 in zip(old_gene_0.split('.'),
                                            old_gene_1.split('.')):
            # Make a random split somewhere inside
            split = random.randint(0, len(gene_part_0))
            new_gene.append(mutate(gene_part_0[:split] + gene_part_1[split:],
                                   mutation_rate))
        if any(x['gene'] == ".".join(new_gene) for x in new_population):
            # Alraedy in the population, don't include twice
            continue
        elif site_replace(structure, replace_list=site_list,
                          manual_angles=new_gene, backends=backends):
            debug("Successfully mated as {}".format(new_gene))
            individual = {'fitness': 0.0,
                          'gene': ".".join(new_gene)}
            new_population.append(individual)
            count -= 1


def mutate(gene, mutation_rate):
    """
    Make a random alteration to a gene string. Randomises letter according to
    the mutation rate.
    """
    new_gene = []
    for allele in gene:
        if random.random() < mutation_rate:
            new_gene.append(random.choice(string.ascii_lowercase))
        else:
            new_gene.append(allele)

    return "".join(new_gene)


def fill_with_randoms(population, max_pop, structure, site_list, backends):
    """
    Make the population up to max_pop in size with completely random
    individuals.

    """
    while len(population) < max_pop:
        this_gene = []
        for functional_group, site in site_list:
            this_gene.append("".join(random.choice(string.ascii_lowercase) for
                                     _ in structure.attachments[site]))
        if site_replace(structure, replace_list=site_list,
                        manual_angles=this_gene, backends=backends):
            individual = {'fitness': 0.0,
                          'gene': ".".join(this_gene)}
            population.append(individual)


def measure_fitness(structure, site_list, individual):
    faps = ['python', '/home/tdaff/bin/faps',
            '-o', 'no_dft', '-o', 'no_force_field_opt=False',
            '-o', 'no_charges', '-o', 'no_gcmc',  '-o', 'no_properties',
            '-o', 'queue=serial', '-o', 'ff_opt_code=gromacs']

    func_name = []
    for site, angles in zip(site_list, individual['gene'].split('.')):
        func_name.append("{}@{}%{}".format(site[0], site[1], angles))

    faps_basename = '%s_func_%s' % (structure.name, ".".join(func_name))
    faps.append(faps_basename)

    faps_run = subprocess.Popen(faps, stdout=subprocess.PIPE)
    faps_run.wait()

    uff_energy = float('inf')
    for line in open("{}.flog".format(faps_basename)):
        if 'UFF energy' in line:
            uff_energy = float(line.split()[-2])

    return uff_energy


def main():
    """
    Main logic to find lowest energy conformation of the functional groups
    through an evolutionary.

    """

    info("Welcome to confswitch; finding your lowest energy conformers")
    info("Using fapswitch version {}".format(fapswitch.__version__))

    # Name for a the single structure
    job_name = options.get('job_name')

    # Load it
    input_structure = load_structure(job_name)
    # Structure is ready!

    # Begin processing
    info("Structure attachment sites: "
         "{}".format(list(input_structure.attachments)))
    info("Structure attachment multiplicities: "
         "{}".format(dict((key, len(val)) for key, val in
                          input_structure.attachments.items())))

    # Functional group library is self initialising
    info("Groups in library: {}".format(functional_groups.group_list))

    # Only use the cif backend here
    backends = [CifFileBackend()]

    # Optimise each custom string passed to fapswitch
    custom_strings = options.get('custom_strings')
    site_strings = re.findall(r'\[(.*?)\]', custom_strings)
    debug("Site replacement options strings: {}".format(site_strings))
    for site_string in site_strings:
        # These should be functional_group1@site1.functional_group2@site2
        # with optional %angle
        site_list = []
        manual_angles = []
        for site in [x for x in site_string.split('.') if x]:
            site_id, functionalisation = site.split('@')
            if '%' in functionalisation:
                functionalisation, manual = functionalisation.split('%')
            else:
                manual = None
            site_list.append([site_id, functionalisation])
            manual_angles.append(manual)

        debug(str(site_list))
        debug(str(manual_angles))

        optimise_conformation(input_structure, site_list, backends=backends)
#                              seed_angles=manual_angles)


if __name__ == '__main__':
    main()
