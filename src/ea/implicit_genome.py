"""
    Data structures and functions associated with the creation
    mapping to neural net (pybrain), mutation, recombination
"""
# -*- coding: utf-8 -*-
#TODO General list: config global params(probas, start and end codons, other)
import time
import random
import sys
import math
import numpy as np
#from pybrain.structure import FeedForwardNetwork, LinearLayer
#from pybrain.structure import FullConnection, SigmoidLayer
import ConfigParser
import itertools
sys.path.insert(0, '../plots')
import plot
import colorama as clr
import exceptions as exc
import brewer2mpl
from pylab import subplot2grid
import load_classif_image as classif
import datetime
import collections
import os


class PrettyFloat(float):
    "Two digit print representation for floating numbers"
    def __repr__(self):
        return "%0.2f" % self


class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


def flatten(l):
    '''From stackoverflow. Transforms a multidimensional list into a flattened
    1-D list'''
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el,
                                                                   basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el


def format_float(value):
    """Four digit print representation for floating numbers"""
    return "%.4f" % value


def map_config(section):
    """Read parameters in given section"""
    dict1 = {}
    options = PARSER.options(section)
    for option in options:
        try:
            dict1[option] = PARSER.get(section, option)
        except exc.IOError:
            print "exception on input/output"
        except Exception:
            print "exception on %s!" % option
            dict1[option] = None
    return dict1


def do_nothing(genome):
    """No fragment mutation"""
    return genome


def random_genome(size, sym):
    """Create random string of given symbols with given size"""
    genome = []
    for _ in xrange(size):
        genome.append(random.choice(sym))
    return genome


def task(vector, task_id, image=None):
    """Return output of given task on given input vector"""
    tasks = {"AND": all, "OR": any, "MAJ": majority, "MIN": minority,
             "RETAND": retina_and, "RETOR": retina_or, "IMG": classify_img}
             # "MUX1": mux1,

    #float vector?
    if task_id == "IMG":
        return classify_img(vector, image)
    else:
        int_vector = [int(round(element)) for element in vector]
        return tasks[task_id](int_vector)


def load_img(filename):
    """External call for loading a png image and returning it as a list of
    pixels starting from top-leftmost pixel, from left to right, then from
    top to bottom"""
    return classif.read_image(filename)


def classify_img(input_x, image):
    """Return the pixel value (0 or 1, for black or white), given coordinates
    in input_x = (x_coord, y_coord)"""
    height = image[1]
    width = image[0]
    y_coordinate = int(round(input_x[1] * height))
    x_coordinate = int(round(input_x[0] * width))
    print y_coordinate, " of ", height, ", and ", x_coordinate, " of ", width

    return list(itertools.chain(image[2]))[y_coordinate][x_coordinate]


def retina_and(input_x):
    """Retina AND problem
    Input vector is arranged in [left1,l2,r1,r2, (2nd row) l3,l4,r3,r4]"""
    l_input = [vector for idx, vector in enumerate(input_x)
               if idx in [0, 1, 4, 5]]
    r_input = [vector for idx, vector in enumerate(input_x)
               if idx in [2, 3, 6, 7]]
    return is_l_object(l_input) and is_r_object(r_input)


def retina_or(input_x):
    """Retina OR problem
    Input vector is arranged in [left1,l2,r1,r2, (2nd row) l3,l4,r3,r4]"""
    l_input = [vector for idx, vector in enumerate(input_x)
               if idx in [0, 1, 4, 5]]
    r_input = [vector for idx, vector in enumerate(input_x)
               if idx in [2, 3, 6, 7]]
    return is_l_object(l_input) or is_r_object(r_input)


def is_l_object(left_x):
    """Is there an object on the left side?"""
    if left_x.count(1) >= 3:
        return True
    else:
        #Only one or two pixeles in left column, none in right column
        if (left_x[1] == 0) and (left_x[3] == 0):
            if (left_x[0] == 1) or (left_x[2] == 1):
                return True
    return False


def is_r_object(right_x):
    """Is there an object on the right side?"""
    if right_x.count(1) >= 3:
        return True
    else:
        #Only one or two pixeles in right column, none in leftcolumn
        if (right_x[0] == 0) and (right_x[2] == 0):
            if (right_x[1] == 1) or (right_x[3] == 1):
                return True
    return False


def majority(input_x):
    """Majority (of half) of ones"""
    total = sum(input_x)
    if total >= len(input_x)/2.0:
        return True
    else:
        return False


def minority(input_x):
    """Minority of ones"""
    return not majority(input_x)
###############################################################################


def frag_copy(genome):
    """Copy a randomly chosen fragment at a randomly chosen position"""
    offspring = []
    size = len(genome)
    #Draw a first position
    first = random.choice(xrange(size + 1))
    #Draw a second position (after "first" position)
    second_range = range(first, size + 1)
    second = random.choice(second_range)
    at_position = random.choice(xrange(size + 1))

    for idx in xrange(size):
        if idx == at_position:
            offspring.extend(genome[first:second])
        offspring.append(genome[idx])

    if at_position == size + 1:
        offspring.extend(genome[first:second])
    return offspring


def frag_move(genome):
    """Move a randomly chosen fragment at a randomly chosen position"""
    offspring = []
    size = len(genome)
    #Draw a first position
    first = random.choice(xrange(size + 1))
    #Draw a second position (after "first" position)
    second_range = range(first, size + 1)
    second = random.choice(second_range)
    frag = genome[first:second]
    #Insert drawn fragment at random position in the pruned offspring
    at_position = random.choice(xrange(len(offspring) + 1))

    offspring = genome[:first] + genome[second:]

    beginning = offspring[0:at_position]
    end = offspring[at_position:]
    offspring = beginning + frag + end

    return offspring


def frag_del(genome):
    """Delete a randomly chosen fragment of the genome"""
    offspring = []
    size = len(genome)
    #Draw interval to delete
    first = random.choice(xrange(size + 1))
    second_range = range(size + 1)[first:size + 1]
    second = random.choice(second_range)
    offspring = genome[:first] + genome[second:]

    return offspring


def bit_flip(genome, prob):
    """Bit-flip at all positions, each with a probability prob
    For binary genomes only.
    Other mutation operators to define in non-binary case"""
    offspring = genome[:]

    for idx, char in enumerate(offspring):
        if np.random.uniform() < prob:
            if char == '0':
                offspring[idx] = '1'
            else:
                if char == '1':
                    offspring[idx] = '0'
                else:
                    sys.exit("Wrong char in genome")
    return offspring


def random_genome_with_genes(junk, nb_links, limit_codons, nucl_per_weight,
                             sym):
    """Create a random genome string with given number of link genes, using
    given start and end codons and given symbol alphabet. It sets a random
    initial weight field of nucl_per_weight characters, and "junk" characters
    in total are added into inter-gene space"""
    start_codon = limit_codons[0]
    end_codon = limit_codons[1]
    genome = []
    interstices = nb_links + 1
    junk_per_interstice = int(math.floor(junk / interstices))
    target_size = max(1, math.ceil(math.log(nb_links, len(sym))))
    #Junk nucleotids at the beginning
    genome = genome + list(np.random.choice(sym, size=junk_per_interstice))
    #loop for each link gene
    for idx in xrange(nb_links):
        #Start codon
        genome = genome + list(start_codon)
        #Targeted link ID in "tgtSize" positions (binary numbering)
        target = list(("{0:0%sb}" % str(int(target_size))).format(idx))

        genome = genome + target
        #Some random nucleotids for weight field
        for _ in xrange(nucl_per_weight):
            genome = genome + list(random.choice(sym))
        #End codon
        genome = genome + list(end_codon)
        #Junk nucleotids before next gene
        genome = genome + list(np.random.choice(sym, size=junk_per_interstice))

    return genome


def create_genome(nb_links, fraction_junk_genes, nucleotids_per_gene,
                  limit_codons, sym):
    """Compute necessary parameters (target field size, effective genome size
    including junk nucleotids), and call random_genome_with_genes to create
    a genome including all specified nb_links"""
    start_codon = limit_codons[0]
    end_codon = limit_codons[1]
    target_size = int(math.ceil(max(1, math.log(nb_links, len(sym)))))
    sze = nb_links * (target_size + len(start_codon) + len(end_codon) +
                      nucleotids_per_gene)
    effective_size = int(round((1.0/(1.0 - fraction_junk_genes)) * sze))

    genome = random_genome_with_genes(
        int(math.floor(effective_size * fraction_junk_genes)),
        nb_links, limit_codons, nucleotids_per_gene, sym)
    #Robustness through redundancy?
    return genome


def map_to_mlp_no_pybrain(codon_list, n_in, n_out, n_hidden=0, neur_per_hid=0,
                          polygene_strategy="avg"):  # bias=False(?)
    """Mapping a codon list (id_codon, weight) to a given neural structure,
    It returns the weight vector(s) connecting each pair of layers in a
    Feed Forward NN, or multilayered perceptron.
    The effect of multiple codons targeting the same
    connection is determined by the polygene_strategy (average by default)"""

    n_weight_layers = n_hidden + 1
    weight_vectors = [[]] * n_weight_layers
    if n_hidden <= 0:
        num_links = n_in * n_out
        weight_vectors[0] = [0.0] * num_links
    else:
        num_links = n_in * neur_per_hid + neur_per_hid * n_out + \
            (n_hidden - 1) * (neur_per_hid * neur_per_hid)
        weight_vectors[0] = [0.0] * n_in * neur_per_hid
        for i in xrange(1, n_weight_layers - 1):
            weight_vectors[i] = [0.0] * neur_per_hid * neur_per_hid
        weight_vectors[-1] = [0.0] * neur_per_hid * n_out

    target_set = set(xrange(num_links))
    grouped_codons = [[codon_id, [codon[1] for codon in codon_list
                                  if codon[0] == codon_id]]
                      for codon_id in target_set]
    grouped_codons = sorted(grouped_codons, key=lambda codon: codon[0])
    if polygene_strategy == "avg":
        weights = dict([[codon[0], sum(codon[1])/float(len(codon[1]))
                         if len(codon[1]) else 0.0]
                        for codon in grouped_codons])
    if n_hidden <= 0:
        for i in xrange(num_links):
            if i in weights:
                weight_vectors[0][i] = weights[i]
        weight_vectors[0] = np.reshape(weight_vectors[0], (n_in, n_out))
    else:
        i = 0
        for l in xrange(n_weight_layers):
            if l == 0:
                for j in xrange(n_in * neur_per_hid):
                    if i in weights:
                        weight_vectors[l][j] = weights[i]
                    i = i + 1
                weight_vectors[l] = np.reshape(weight_vectors[l],
                                               (n_in, neur_per_hid))
            else:
                if l == (n_weight_layers - 1):
                    for j in xrange(neur_per_hid * n_out):
                        if i in weights:
                            weight_vectors[l][j] = weights[i]
                        i = i + 1
                    weight_vectors[l] = np.reshape(weight_vectors[l],
                                                   (neur_per_hid, n_out))
                else:
                    for j in xrange(neur_per_hid * neur_per_hid):
                        if i in weights:
                            weight_vectors[l][j] = weights[i]
                        i = i + 1
                    weight_vectors[l] = np.reshape(weight_vectors[l],
                                                   (neur_per_hid,
                                                    neur_per_hid))
    return weight_vectors


def xtract_codons(genome, target_size, limit_codons, max_weight=1.0,
                  pleio=False):
    """Extracts the sequence of coding genes from genome, using given start
    and end codons and with given targeted link string size. Connection weight
    is capped in [-max_weight, max_weight]. Pleiotropy can be switched on and
    off. It also returns a user-friendly representation of the genomes
    with segmentated and color-coded codons """
    start_codon = limit_codons[0]
    end_codon = limit_codons[1]
    idx = 0
    out_str = ""
    gene_str = ""
    spacing = "  ||  "
    separator = " - "
    coding_part_size = 0
    codon_limit_list = []
    while idx in xrange(len(genome)):
        #If there is a start codon, then read a whole gene
        if "".join(genome[idx:]).startswith(start_codon):
            gene_start_index = idx
            idx = idx + len(start_codon)

            gene_str = spacing + clr.Back.YELLOW \
                + "".join(genome[gene_start_index:idx]) +\
                clr.Back.RESET + separator
            #Read target field: it corresponds to the targeted link ID
            tgt = genome[idx:idx + target_size]
            gene_str = gene_str + clr.Back.CYAN + "".join(tgt) +\
                clr.Back.RESET + separator
            idx = idx + target_size
            j = idx
            #Look for the end codon, in order to delimitate the weight field:
            #the field storing the weight
            while not("".join(genome[j:]).startswith(end_codon)) and\
                    j < len(genome):
                j = j + 1
            #If there was an end codon, then the gene was valid, so we add it
            #to the list of codons (ID, weight)
            if j < len(genome):
                weight_string = genome[idx:j]
                gene_str = gene_str + clr.Back.GREEN +\
                    "".join(weight_string) + clr.Back.RESET + separator

                gene_str = gene_str + clr.Back.YELLOW + \
                    "".join(genome[j:j + len(end_codon)]) +\
                    clr.Back.RESET + spacing

                out_str = out_str + gene_str
                gene_str = ""

                codon_limit_list.append([gene_start_index, j +
                                         len(end_codon)])
            else:
                #if genome read is finished and the current gene is not ended
                #it is not a valid gene, and it is not added to the codon list
                out_str = out_str + "".join(genome[gene_start_index:])
                gene_str = ""
                break
            if pleio:
                #If pleiotropy is activated, continue reading just after
                #previous gene's start codon
                idx = idx - target_size
            else:
                #Continue reading after end codon otherwise
                idx = j + len(end_codon)
        #If there is not a start codon, keep reading
        else:
            out_str = out_str + str(genome[idx])
            idx = idx + 1
    #Measuring coding part
    coding_part_size = 0
    sum_per_gene_coding_size = 0

    if len(codon_limit_list) > 0:
        i = 0
        start_gene = codon_limit_list[i][0]
        end_gene = codon_limit_list[i][1]
        sum_per_gene_coding_size = sum_per_gene_coding_size +\
            codon_limit_list[i][1] - codon_limit_list[i][0]
        i = i + 1
        finished = (i == len(codon_limit_list))
        while not finished:
            codon_limits = codon_limit_list[i]
            l1 = codon_limits[0]
            l2 = codon_limits[1]
            if l1 >= end_gene:
                #separated gene
                coding_part_size = coding_part_size + end_gene - start_gene
                start_gene = l1
                end_gene = l2
            else:
                end_gene = l2
            i = i + 1
            finished = (i == len(codon_limit_list))
            sum_per_gene_coding_size = sum_per_gene_coding_size + l2 - l1
        #add last gene's size
        coding_part_size = coding_part_size + end_gene - start_gene

    codon_list = []

    for codon_limits in codon_limit_list:
        gene = genome[codon_limits[0]:codon_limits[1]]
        tgt = gene[len(start_codon):len(start_codon) + target_size]
        weight_string = gene[len(start_codon) + target_size:-len(end_codon)]
        if len(weight_string) == 0:
            codon_list.append([int("".join(tgt), 2), 0.0])
        else:
            codon_list.append([int("".join(tgt), 2), 2.0 * max_weight *
                               weight_string.count("1") /
                               len(weight_string) - max_weight])

    return codon_list, out_str, coding_part_size, sum_per_gene_coding_size
###############################################################################


def mutate(genome, probs):
    """Mutate one single genome.
    First, copy, move or delete a random fragment of the genome at random pos.
    Second, bit-flip for each bit of the genome with small probability"""

    p_frag_copy = probs[0]
    p_frag_move = probs[1]
    p_frag_del = probs[2]
    p_bit_flip = probs[3]

    muts = {0: frag_copy,
            1: frag_move,
            2: frag_del,
            3: do_nothing}

    p_nothing = 1 - p_frag_copy - p_frag_move - p_frag_del

    mut = np.random.choice(xrange(0, 4),
                           p=[p_frag_copy, p_frag_move, p_frag_del, p_nothing])

    offspring = bit_flip(muts[mut](genome), p_bit_flip)

    return offspring


def stats_codons(l_codons, genome, n_links):
    '''Computes some statistics on an individual, mainly regarding coding genes
       in the genome, cf. below for details.
       Returns: (number of coding genes, average number of genes for the same
       link, size of the coding part in the genotype)'''
    result = [0.0, 0.0, 0.0]
    target_set = set(xrange(n_links))
    grouped_codons = [[codon_id, [codon[1] for codon in l_codons
                                  if codon[0] == codon_id]]
                      for codon_id in target_set]
    codons_per_link = [len(link[1]) for link in grouped_codons]

    result[0] = float(len(l_codons))
    result[1] = np.average(codons_per_link)
    result[2] = float(len(genome))
    return result


def diversity(pop, nb_link):
    '''Computes measures of diversity in the population based on the
       differences between all pairs of individuals.'''

    vector_stats = []
    centroid = [0.0] * nb_link
    for indiv in pop:
        centroid = np.add(centroid, list(flatten(indiv[4])))
    for i in xrange(len(centroid)):
        centroid[i] = centroid[i] / len(pop)

    for i, indiv in enumerate(pop):
        stats = stats_codons(indiv[2], indiv[0], nb_link)
        target_set = set(xrange(nb_link))
        grouped_codons = [[codon_id, [codon[1] for codon in indiv[2]
                                      if codon[0] == codon_id]]
                          for codon_id in target_set]
        nb_valid_links = len(grouped_codons)
        if nb_valid_links != 0:
            weight_norm = np.linalg.norm(list(flatten(indiv[4]))) /\
                nb_valid_links
        else:
            weight_norm = 0.0

        if len(indiv[0]) != 0:
            coding_frac = float(indiv[5]) / float(len(indiv[0]))
        else:
            coding_frac = 0.0

        if float(indiv[6]) != 0.0:
            deg_pleio = 1 - (float(indiv[5]) / float(indiv[6]))
        else:
            deg_pleio = 0.0
        if deg_pleio > 1.0 or deg_pleio < 0.0:
            print "Pleio. deg: ", deg_pleio
        vector_stats.append(stats +
                            [coding_frac, float(indiv[5]), weight_norm,
                             np.linalg.norm(np.subtract(list(flatten(indiv[4])),
                                            centroid)), deg_pleio])
    return vector_stats
    #number of codons, avg number of genes per link, genome size, coding genome
    #fraction, coding genome size, weight vector norm


def evaluate(ind, problem_db, lbl_db):
    """Evaluate indiv NN on nb_instances
    with input=problem_db and output=lbl_db"""
    fit = 0

    for idx, inp in enumerate(problem_db):
        nn_answer = ind.activate(inp)
        #Error on each instance (euclidian distance to right answer)
        #is normalized between 0 and 1
        error = np.linalg.norm(np.subtract(nn_answer, lbl_db[idx])) / \
            math.sqrt(len(inp))
        fit += 1.0 - error
    #Normalized for the size of the database:
    #average normalized euclidian error over instances
    return fit/float(len(problem_db))


def activate(vector_weights, in_pattern):
    '''Computes and returns the output activation of a layered FF neural
    network having the given vector of weights and the input pattern'''

    n_layers = len(vector_weights)
    activation = in_pattern

    for i in xrange(n_layers):
        activation = np.dot(np.array(activation), vector_weights[i])
        activation = np.tanh(activation)
    return np.array(activation)


def evaluate_no_pybrain(ind, problem_db, lbl_db):
    """Evaluate indiv NN on nb_instances
    with input=problem_db and output=lbl_db"""
    fit = 0

    for idx, inp in enumerate(problem_db):
        nn_answer = activate(ind, inp)
        #Error on each instance (euclidian distance to right answer)
        #is normalized between 0 and 1
        error = np.linalg.norm(np.subtract(nn_answer, lbl_db[idx])) / \
            math.sqrt(len(inp))
        fit += 1.0 - error
    #Normalized for the size of the database:
    #average normalized euclidian error over instances
    return fit/float(len(problem_db))


def select(population, nb_parents):
    """Parent rank-based selection of lambda (= nb_parents) parents, to be
    mutated afterwards"""
    selected_indexes = [("", 0.0, [], "", None)] * nb_parents
    fitness_list = [ind[1] for ind in population]
    sorted_indexes = [idx[0] for idx in sorted(enumerate(fitness_list),
                                               key=lambda x: x[1])]
    rank_weight = [idx[0] for idx in enumerate(sorted_indexes)]
    triang_number = len(fitness_list) * (len(fitness_list) + 1)/2
    rank_proba = [float(idx + 1)/float(triang_number) for idx in rank_weight]
    #Select with replacement with rank-based probabilities
    for idx in xrange(nb_parents):
        chosen = np.random.choice(sorted_indexes, p=rank_proba)
        selected_indexes[idx] = chosen

    condition = (nb_parents == len(selected_indexes))
    assert condition, "Wrong number of selected parents"
    selected_parents = [population[ind] for ind in selected_indexes]
    return selected_parents


def survive(population, mu_nb_par):
    """Truncation + survivor selection of the best mu_nb_par individuals
    of both children and parents"""
    sorted_individuals = [ind for ind in sorted(population,
                                                key=lambda x: x[1])]
    survivors = sorted_individuals[-mu_nb_par:]
    return survivors


if __name__ == "__main__":
    clr.init()
    CONFIG = sys.argv[1]
    OUT_FOLDER = sys.argv[2]
    PARSER = ConfigParser.ConfigParser()
    PARSER.read(CONFIG)
    PARAMS = AutoVivification()

    for s in PARSER.sections():
        PARAMS[s] = map_config(s)

    VERBOSE = bool(PARAMS['EA']['verbose'])

    TASK_SEQUENCE_LIST = PARAMS["Task"]["tasksequence"].split(",")
    OUT_FOLDER_NAME = CONFIG.split("/")[-1].split(".")[0]

    if not os.path.exists("logs/" + OUT_FOLDER_NAME):
        os.mkdir("logs/" + OUT_FOLDER_NAME)

    if not os.path.exists("logs/" + OUT_FOLDER_NAME + "/" + OUT_FOLDER):
        os.mkdir("logs/" + OUT_FOLDER_NAME + "/" + OUT_FOLDER)

    LOG_FILENAME = "logs/" + OUT_FOLDER_NAME + "/" +\
        OUT_FOLDER + "/fitness_" +\
        str(datetime.datetime.now()).replace(" ", "-") + ".log"
    LOG_FILE = open(LOG_FILENAME, 'w')
    TEST_LOG_FILENAME = "logs/" + OUT_FOLDER_NAME + "/" +\
        OUT_FOLDER + "/test_" +\
        str(datetime.datetime.now()).replace(" ", "-") + ".log"
    TEST_LOG_FILE = open(TEST_LOG_FILENAME, 'w')
    NN_LOG_FILENAME = "logs/" + OUT_FOLDER_NAME + "/" +\
        OUT_FOLDER + "/neural_" +\
        str(datetime.datetime.now()).replace(" ", "-") + ".log"
    NN_LOG_FILE = open(NN_LOG_FILENAME, 'w')
    SIZE_LOG_FILENAME = "logs/" + OUT_FOLDER_NAME + "/" +\
        OUT_FOLDER + "/size_" +\
        str(datetime.datetime.now()).replace(" ", "-") + ".log"
    SIZE_LOG_FILE = open(SIZE_LOG_FILENAME, 'w')
    CODONS_LOG_FILENAME = "logs/" + OUT_FOLDER_NAME +\
        "/" + OUT_FOLDER + "/codons_" +\
        str(datetime.datetime.now()).replace(" ", "-") + ".log"
    CODONS_LOG_FILE = open(CODONS_LOG_FILENAME, 'w')
    GENES_PER_LINK_LOG_FILENAME = "logs/" +\
        OUT_FOLDER_NAME +\
        "/" + OUT_FOLDER + "/gpl_" +\
        str(datetime.datetime.now()).replace(" ", "-") + ".log"
    GENES_PER_LINK_LOG_FILE = open(GENES_PER_LINK_LOG_FILENAME, 'w')
    CODING_SIZE_LOG_FILENAME = "logs/" +\
        OUT_FOLDER_NAME +\
        "/" + OUT_FOLDER + "/coding_size_" +\
        str(datetime.datetime.now()).replace(" ", "-") + ".log"
    CODING_SIZE_LOG_FILE = open(CODING_SIZE_LOG_FILENAME, 'w')
    WEIGHT_NORM_LOG_FILENAME = "logs/" +\
        OUT_FOLDER_NAME +\
        "/" + OUT_FOLDER + "/weight_norm_" +\
        str(datetime.datetime.now()).replace(" ", "-") + ".log"
    WEIGHT_NORM_LOG_FILE = open(WEIGHT_NORM_LOG_FILENAME, 'w')
    PLEIO_LOG_FILENAME = "logs/" +\
        OUT_FOLDER_NAME +\
        "/" + OUT_FOLDER + "/pleio_" +\
        str(datetime.datetime.now()).replace(" ", "-") + ".log"
    PLEIO_LOG_FILE = open(PLEIO_LOG_FILENAME, 'w')
###############################################################################
    # 1 output for logical binary problems
    N_OUT = int(PARAMS["Task"]["outputs"])
    N_IN = int(PARAMS["Task"]["inputs"])  # 2  # 8 for retina pb
    N_HID = int(PARAMS['NN']['nhidlay'])  # 1
    NEUR_HID_LAYER = int(PARAMS['NN']['neurperlay'])  # 2
    #Proportion of junk genes in-between genes on initialization of the genome
    FRAC_JUNK_GENES = float(PARAMS['Encoding']['junk'])
    SYMBOLS = ['0', '1']
    #Number of nucleotids for weight field on initialization
    #(the more nucleotids, the finer the resolution of initial weights)
    NB_WEIGHT_VALUES = int(PARAMS['Encoding']['nweightsvalues'])
    NUCL_PER_WEIGHT = NB_WEIGHT_VALUES - 1
    #How to select start and end codons? Shorter (i.e. easier to randomly draw)
    #codons are more prone to disruptive variation? (are offspring viable?)
    #Grey-coding for links? minimizing "distances" between links?
    #Look for a symmetry in codons?

    START_CODON = PARAMS['Encoding']['start']
    END_CODON = PARAMS['Encoding']['end']

    PLEIOTROPY = bool(PARAMS['Encoding']['pleio'])
    if N_HID == 0:
        N_LINKS = N_IN * N_OUT
    else:
        N_LINKS = N_IN * NEUR_HID_LAYER + (N_HID - 1) * NEUR_HID_LAYER**2 +\
            NEUR_HID_LAYER * N_OUT
    TGT_SIZE = int(math.ceil(max(1, math.log(N_LINKS, len(SYMBOLS)))))
    hid_shape = [NEUR_HID_LAYER for _ in xrange(N_HID)]
    nn_shape = tuple([N_IN] + hid_shape + [N_OUT])

#############################Evolutionary parameters###########################
    MU = int(PARAMS['EA']['mu'])
    NB_GENERATIONS = int(PARAMS['EA']['generations'])
    # number of children per parent
    PARENT_CHILDREN_RATIO = float(PARAMS['EA']['lambdapermu'])
    # parent/children ratio * mu == lambda nb of children
    LAMBDA = int(round(MU * PARENT_CHILDREN_RATIO))
    #mu parents in the beginning and lbda children every generation
    mut_prob = []
    mut_prob.append(float(PARAMS['MutProb']['pcpy']))
    mut_prob.append(float(PARAMS['MutProb']['pmove']))
    mut_prob.append(float(PARAMS['MutProb']['pdel']))
    mut_prob.append(float(PARAMS['MutProb']['pflip']))

###############################################################################
    PROBLEM_DB = []
    LABEL_DB = []
    PROBLEM_SIZE = N_IN
    #-1  # -1 for whole problem
    TASK_SEQUENCE = TASK_SEQUENCE_LIST

    PROBLEM_ID = PARAMS['Task']['problemid']
    # "RETAND" "MIN", "AND", "OR", "MAJ", "RETOR" "IMG"
###############################################################################
###############################################################################
##########################Evolutionary algorithm###############################
    INTERTASK_LOG = []
    POPULATION = [("", 0.0, [], "", None, 0.0, 0.0, 0.0)] * MU

    #Initial genomes in population: mu individuals, valid ones
    #(all random links correctly encoded)
    for i in xrange(MU):
        g = create_genome(N_LINKS, FRAC_JUNK_GENES, NUCL_PER_WEIGHT,
                          (START_CODON, END_CODON), SYMBOLS)
        codons, PRETTY_STRING, c_size, sum_genes_size =\
            xtract_codons(g, TGT_SIZE, (START_CODON, END_CODON),
                          pleio=PLEIOTROPY)
        #Only consider codons with a valid (existing) targeted link
        codons = [c for c in codons if c[0] < N_LINKS]
        net = map_to_mlp_no_pybrain(codons, N_IN, N_OUT, n_hidden=N_HID,
                                    neur_per_hid=NEUR_HID_LAYER)
        individual = (g, 0.0, codons, PRETTY_STRING, net, c_size,
                      sum_genes_size, 0.0)
        POPULATION[i] = individual

###############################################################################
    #Total Time
    TIME_START = time.time()
    for each_task in TASK_SEQUENCE:
        POPULATION_LOG = []
        DIV_LOG = []

        print "Learning task: ", each_task
        #Generate full problem
        if PROBLEM_ID == "IMG":
            IMG_FILENAME = each_task
            IMG = load_img(IMG_FILENAME)
            PROBLEM_DB = [[y / float(IMG[1]), x / float(IMG[0])]
                          for y in xrange(IMG[1]) for x in xrange(IMG[0])]
            LABEL_DB = [x / 255 for x in IMG[2]]
            TOTAL_INSTANCES = IMG[0] * IMG[1]
        else:
            #Generate all instances
            TOTAL_INSTANCES = len(SYMBOLS)**PROBLEM_SIZE
            for i in xrange(TOTAL_INSTANCES):
                #binary input patterns
                input_vector = [float(k) for k in
                                list(("{0:0%sb}" %
                                      str(PROBLEM_SIZE)).format(i))]
                PROBLEM_DB.append(input_vector)
                LABEL_DB.append(task(input_vector, each_task))

        INSTANCES_DB = []
        LABEL_INST_DB = []
        INSTANCES_TEST_DB = []
        LABEL_TEST_DB = []
        DB_SIZE = int(float(PARAMS['Task']['trainingfraction']) *
                      TOTAL_INSTANCES)

        if DB_SIZE != -1:
            index_instances = np.random.choice(xrange(TOTAL_INSTANCES),
                                               DB_SIZE, replace=False)
            index_test = np.random.choice(xrange(TOTAL_INSTANCES),
                                          DB_SIZE * 5, replace=False)
                                        # TODO how many test examples (5 x sz?)
            INSTANCES_DB = [PROBLEM_DB[i] for i in index_instances]
            LABEL_INST_DB = [LABEL_DB[i] for i in index_instances]
            INSTANCES_TEST_DB = [PROBLEM_DB[i] for i in index_test]
            LABEL_TEST_DB = [LABEL_DB[i] for i in index_test]
        else:
            INSTANCES_DB = PROBLEM_DB[:]
            LABEL_INST_DB = LABEL_DB[:]
            INSTANCES_TEST_DB = PROBLEM_DB[:]
            LABEL_TEST_DB = LABEL_DB[:]
        TIME_GEN_LOG = []

        PRETTY_STRING = ""
        #Initial population: mu individuals,
        #all random links correctly encoded?
        for i in xrange(MU):
            net = POPULATION[i][4]
            #Evaluate initial population
            fitness = evaluate_no_pybrain(net, INSTANCES_DB, LABEL_INST_DB)
            test_perf = evaluate_no_pybrain(net, INSTANCES_TEST_DB,
                                            LABEL_TEST_DB)
            POPULATION[i] = (POPULATION[i][0], fitness, POPULATION[i][2],
                             POPULATION[i][3], POPULATION[i][4],
                             POPULATION[i][5], POPULATION[i][6], test_perf)

        POPULATION_LOG.append(POPULATION)
        DIV_LOG.append(diversity(POPULATION, N_LINKS))
        ITERATION = 0
        sys.stdout.write(str(ITERATION))

        #Evolutionary loop
        for ITERATION in xrange(1, NB_GENERATIONS + 1):
            #Time per generation
            time_gen_start = time.time()
            #Select lambda (possibly repeated) parents, rank-based
            parents = select(POPULATION, LAMBDA)

            children = [("", 0.0, [], "", None, 0.0, 0.0, 0.0)] * LAMBDA

            #Generate lambda children by mutating lambda selected parents
            #Look for disruptive mutations (e.g. by measuring valid codons)
            #Note: it is easier to break a link gene than creating a new one
            for index, i in enumerate(parents):
                child = mutate(i[0][:], mut_prob)
                codons, PRETTY_STRING, c_size, sum_genes_size =\
                    xtract_codons(child[:], TGT_SIZE,
                                  (START_CODON, END_CODON),
                                  pleio=PLEIOTROPY)
                #Only consider codons with a valid (existing) targeted link
                codons = [c for c in codons if c[0] < N_LINKS]

                net = map_to_mlp_no_pybrain(codons, N_IN, N_OUT,
                                            n_hidden=N_HID,
                                            neur_per_hid=NEUR_HID_LAYER)
                #Evaluate children
                fitness = evaluate_no_pybrain(net, INSTANCES_DB, LABEL_INST_DB)
                test_perf = evaluate_no_pybrain(net, INSTANCES_TEST_DB,
                                                LABEL_TEST_DB)
                individual = (child[:], fitness, codons, PRETTY_STRING,
                              net, c_size, sum_genes_size, test_perf)
                children[index] = individual

            full_population = (POPULATION + children)
            #Truncation "plus" survivor sel./replacement [parents + children]
            surviving_individuals = survive(full_population, MU)

            POPULATION = [("", 0.0, [], "", None, 0.0, 0.0, 0.0)] * MU
            for index, individual in enumerate(surviving_individuals):
                POPULATION[index] = individual

            POPULATION_LOG.append(POPULATION)
            DIV_LOG.append(diversity(POPULATION, N_LINKS))

            #Logging operations at the end of generation
            sys.stdout.write(" - " + str(ITERATION))
            time_gen_end = time.time()
            TIME_GEN_LOG.append([time_gen_end - time_gen_start])

        INTERTASK_LOG.append(POPULATION_LOG)
        formatted_fitness = [[format_float(indiv[1]) for indiv in generation]
                             for generation in POPULATION_LOG]
        formatted_test_perf = [[format_float(indiv[7]) for indiv in generation]
                               for generation in POPULATION_LOG]
        formatted_size = [[format_float(indiv[2]) for indiv in generation]
                          for generation in DIV_LOG]
        formatted_codons = [[format_float(indiv[0]) for indiv in generation]
                            for generation in DIV_LOG]
        formatted_genes_per_link = [[format_float(indiv[1])
                                     for indiv in generation]
                                    for generation in DIV_LOG]
        formatted_coding_size = [[format_float(indiv[3])
                                 for indiv in generation]
                                 for generation in DIV_LOG]
        formatted_weight_norm = [[format_float(indiv[6])  # 5 for weight norm
                                  for indiv in generation]
                                 for generation in DIV_LOG]
        formatted_pleio = [[format_float(indiv[7])
                            for indiv in generation]
                           for generation in DIV_LOG]
        LOG_FILE.write(str(formatted_fitness) + "\n")
        TEST_LOG_FILE.write(str(formatted_test_perf) + "\n")
        SIZE_LOG_FILE.write(str(formatted_size) + "\n")
        CODONS_LOG_FILE.write(str(formatted_codons) + "\n")
        GENES_PER_LINK_LOG_FILE.write(str(formatted_genes_per_link) + "\n")
        CODING_SIZE_LOG_FILE.write(str(formatted_coding_size) + "\n")
        WEIGHT_NORM_LOG_FILE.write(str(formatted_weight_norm) + "\n")
        PLEIO_LOG_FILE.write(str(formatted_pleio) + "\n")

        #transfer to next task: best individual of last generation
        max_fitness = POPULATION[0][1]
        best_individual = POPULATION[0]
        for individual in POPULATION:
            if individual[1] > max_fitness:
                max_fitness = individual[1]
                best_individual = individual
        NN_LOG_FILE.write(each_task + "\n" + "".join(best_individual[0]) +
                          "\n" + str(format_float(best_individual[1])) + "\n")
        #If there is a following task, i.e. current is not last task
        if each_task != TASK_SEQUENCE[-1]:
            #Initialize population with mutated copies of best_individual
            for index in xrange(MU):
                altered_copy = mutate(best_individual[0][:], mut_prob)
                codons, PRETTY_STRING, c_size, sum_genes_size =\
                    xtract_codons(altered_copy, TGT_SIZE,
                                  (START_CODON, END_CODON), pleio=PLEIOTROPY)
                #Only consider codons with a valid (existing) targeted link
                codons = [c for c in codons if c[0] < N_LINKS]
                net = map_to_mlp_no_pybrain(codons, N_IN, N_OUT,
                                            n_hidden=N_HID,
                                            neur_per_hid=NEUR_HID_LAYER)
                individual = (altered_copy, 0.0, codons,
                              PRETTY_STRING, net, c_size, sum_genes_size, 0.0)
                POPULATION[index] = individual
        print
###############################################################################
    LOG_FILE.close()
    TEST_LOG_FILE.close()
    NN_LOG_FILE.close()
    SIZE_LOG_FILE.close()
    CODONS_LOG_FILE.close()
    GENES_PER_LINK_LOG_FILE.close()
    CODING_SIZE_LOG_FILE.close()
    WEIGHT_NORM_LOG_FILE.close()
    PLEIO_LOG_FILE.close()
    if VERBOSE:
        TIME_END = time.time()
        print "\nIt took: ", str(TIME_END - TIME_START), " seconds"
        print "Time per generation (seconds)"
        plot.draw_data([["Time", [list(x) for x in zip(*TIME_GEN_LOG)]]])
        print "For mu = ", str(MU), ", lambda = ", str(LAMBDA), \
            ", number of generations = ", NB_GENERATIONS, ", problem = ",\
            PROBLEM_ID
        print "\nFitness through generations (1 - error)"
        FITNESS_LOG = [[f for _, f, _, _, _ in generation]
                       for generation in POPULATION_LOG]
        plot.draw_data([["Test", [list(x) for x in zip(*FITNESS_LOG)]]])
        print "Number of codons (of the mu survirving individuals) \
                                                                per generation"
        NB_CODONS_LOG = [[len(indiv[2]) for indiv in generation]
                         for generation in POPULATION_LOG]
        plot.draw_data([["Codons",
                         [list(x) for x in zip(*NB_CODONS_LOG)]]])

###############################################################################
        #Test (& generalization if tested on more instances)
        NETS_POPULATION = []
        for g in POPULATION:
            codons, PRETTY_STRING, c_size, sum_genes_size =\
                xtract_codons(g[0], TGT_SIZE, (START_CODON, END_CODON),
                              pleio=PLEIOTROPY)
            #Only consider codons with a valid (existing) targeted link
            codons = [c for c in codons if c[0] < N_LINKS]
            net = map_to_mlp_no_pybrain(codons, N_IN, N_OUT, n_hidden=N_HID,
                                        neur_per_hid=NEUR_HID_LAYER)
            NETS_POPULATION.append(net)
        TEST_RESULTS = []
        for net in NETS_POPULATION:
            accuracy = 0.0
            for i in xrange(2**PROBLEM_SIZE):
                answer = net.activate(PROBLEM_DB[i])
                if int(round(answer)) == LABEL_DB[i]:
                    accuracy = accuracy + 1.0
            accuracy = accuracy / (2 ** PROBLEM_SIZE)
            TEST_RESULTS.append(accuracy)
        print "Test results (accuracy of last generation)"
        AXIS = subplot2grid((1, 1), (0, 0))
        BMAP = brewer2mpl.get_map('Set2', 'qualitative', 7)
        plot.plot_boxplot([TEST_RESULTS], BMAP.mpl_colors, AXIS, ["Test"])
