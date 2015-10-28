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
from pybrain.structure import FeedForwardNetwork, LinearLayer
from pybrain.structure import FullConnection, SigmoidLayer
import ConfigParser
import itertools
sys.path.insert(0, '../plots')
import plot
import colorama as clr
import exceptions as exc
import brewer2mpl
from pylab import subplot2grid
import load_classif_image as classif
#import parse; import difflib
import datetime


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


def format_float(value):
    """Three digit print representation for floating numbers"""
    return "%.3f" % value


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


def map_to_standard_mlp(codon_list, n_in, n_out, n_hidden=0, neur_per_hid=0,
                        polygene_strategy="avg"):  # bias=False(?)
    """Map a list of coding genes (codons) to a Multilayer Fully-connected
    Perceptron with given inputs, outputs, number of hidden layers and neurons
    per hidden layer. The effect of multiple codons targeting the same
    connection is determined by the polygene_strategy (average by default)"""
    mapped_net = None
    if n_hidden <= 0:
        num_links = n_in * n_out
    else:
        num_links = n_in * neur_per_hid + neur_per_hid * n_out + \
            (n_hidden - 1) * (neur_per_hid * neur_per_hid)
    #Neurons' IDs starts counting from 0 in the sequence:
    #Input - Bias - Output - (Hidden)
    #Neuron IDs are fixed and topology does not evolve for now
    #Link order : cf. below
    #Recurrent NNs? evolution of topology?

    mapped_net = FeedForwardNetwork()

    for idx in xrange(n_in):
        mapped_net.addInputModule(LinearLayer(1, name="i" + str(idx)))

    for idx in xrange(n_out):
        mapped_net.addOutputModule(SigmoidLayer(1, name="o" + str(idx)))

    for layer in xrange(n_hidden):
        for idx in xrange(neur_per_hid):
            mapped_net.addModule(SigmoidLayer(1, name="h" + str(layer)
                                              + "-" + str(idx)))

    #[[Bias??]] if bias: b = BiasUnit(name = "b");  mapped_net.addModule(b)

    #Loop link codons
    #Link ordering:
    #    No hidden layer:
    #              0        1     ....
    #           [i0,o0], [i1,o0], .... , [i0,o1], [i1,o1], ...., [iN,oN']
    #   W. hidden layer(s):
    #           [i0,h00],[i1,h00],...,[i0,h01],...
    #           [h00,h10],[h01,h10],..., [hK0,o0], [hK1,o0],..., [hKN'',oN']
    input_idx = 0  # for the input neuron
    output_idx = 0  # for the output neuron
    h_layer = 0  # for looping over hidden layers

    target_set = set(xrange(num_links))
    grouped_codons = [[codon_id, [codon[1] for codon in codon_list
                                  if codon[0] == codon_id]]
                      for codon_id in target_set]
    #Loop over all possible link identifiers
    for link in xrange(num_links):
        #Retrieve the list of weights of codons pointing to current link
        gene_weights = [codon[1] for codon
                        in grouped_codons if codon[0] == link]
        gene_weights = list(itertools.chain.from_iterable(gene_weights))

        #if there is no codon pointing to current link, the link is still
        #created, but the weight is 0.0
        if len(gene_weights) == 0:
            weight = 0.0
        else:
            #Different polygenic schemes, dominance etc...
            #i.e. how to integrate several weights of codons pointing to the
            #same link into a single weight to assign to the link of neural net
            #Average
            if polygene_strategy == "avg":
                weight = sum(gene_weights)/float(len(gene_weights))
            #other schemes?
        #Compute the input and output neurons for current link
        #Single Perceptron w/o hidden layer
        if n_hidden == 0:
            in_neuron = mapped_net["i" + str(input_idx)]
            out_neuron = mapped_net["o" + str(output_idx)]
            if input_idx == n_in - 1:
                #To first in neuron for next out neuron
                input_idx = 0
                if output_idx == n_out - 1:
                    #Finished setting up all links (perceptron w/o hidd. layer)
                    output_idx = 0
                    conn = FullConnection(in_neuron, out_neuron,
                                          name=str(link))
                    mapped_net.addConnection(conn)
                    conn.params[0] = weight
                    break
                else:
                    #Next out neuron
                    output_idx = output_idx + 1
            else:
                #Next in neuron
                input_idx = input_idx + 1
        #Multi-layer perceptron
        else:
            if h_layer == 0:
                in_neuron = mapped_net["i" + str(input_idx)]
                out_neuron = mapped_net["h" + str(h_layer) +
                                        "-" + str(output_idx)]
                if input_idx == n_in - 1:
                    input_idx = 0
                    if output_idx == neur_per_hid - 1:
                        output_idx = 0
                        h_layer = h_layer + 1
                    else:
                        output_idx = output_idx + 1
                else:
                    input_idx = input_idx + 1
            else:
                if h_layer < n_hidden:
                    in_neuron = mapped_net["h" + str(h_layer - 1) + "-" +
                                           str(input_idx)]
                    out_neuron = mapped_net["h" + str(h_layer)
                                            + "-" + str(output_idx)]
                    if input_idx == neur_per_hid - 1:
                        input_idx = 0
                        if output_idx == neur_per_hid - 1:
                            output_idx = 0
                            h_layer = h_layer + 1
                        else:
                            output_idx = output_idx + 1
                    else:
                        input_idx = input_idx + 1
                else:
                    in_neuron = mapped_net["h" + str(h_layer - 1) + "-" +
                                           str(input_idx)]
                    out_neuron = mapped_net["o" + str(output_idx)]
                    if input_idx == neur_per_hid - 1:
                        input_idx = 0
                        if output_idx == n_out - 1:
                            output_idx = 0
                            conn = FullConnection(in_neuron, out_neuron,
                                                  name=str(link))
                            mapped_net.addConnection(conn)
                            conn.params[0] = weight
                            break
                        else:
                            output_idx = output_idx + 1
                    else:
                        input_idx = input_idx + 1
        conn = FullConnection(in_neuron, out_neuron, name=str(link))
        mapped_net.addConnection(conn)
        conn.params[0] = weight

    mapped_net.sortModules()
    return mapped_net


def extract_codons(genome, target_size, limit_codons, max_weight=1.0,
                   pleio=False):
    """Extracts the sequence of coding genes from genome, using given start
    and end codons and with given targeted link string size. Connection weight
    is capped in [-max_weight, max_weight]. Pleiotropy can be switched on and
    off. It also returns a user-friendly representation of the genomes
    with the codons segmentated and color-coded"""
    start_codon = limit_codons[0]
    end_codon = limit_codons[1]
    codon_list = []
    idx = 0
    out_str = ""
    gene_str = ""
    spacing = "  ||  "  # " \n---\n "
    separator = " - "

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
                if len(weight_string) == 0:
                    codon_list.append([int("".join(tgt), 2), 0.0])
                else:
                    codon_list.append([
                        int("".join(tgt), 2),
                        2.0 * max_weight * weight_string.count("1") /
                        len(weight_string) - max_weight])

                gene_str = gene_str + clr.Back.YELLOW + \
                    "".join(genome[j:j + len(end_codon)]) +\
                    clr.Back.RESET + spacing

                out_str = out_str + gene_str
                gene_str = ""
            else:
                #if genome read is finished and the current gene is not ended
                #it is not a valid gene, and it is not added to the codon list
                out_str = out_str + "".join(genome[gene_start_index:])
                gene_str = ""
                break
            #TOTEST pleiotropy
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
    return codon_list, out_str
###############################################################################


def mutate(genome):
    """Mutate one single genome.
    First, copy, move or delete a random fragment of the genome at random pos.
    Second, bit-flip for each bit of the genome with small probability"""
    #0.12; 0.12; 0.12; 0.01 MIN
    p_frag_copy = 0.05
    p_frag_move = 0.05
    p_frag_del = 0.05
    p_bit_flip = 0.01

    muts = {0: frag_copy,
            1: frag_move,
            2: frag_del,
            3: do_nothing}  # 4 : [geneIns, pGeneIns]};

    p_nothing = 1 - p_frag_copy - p_frag_move - p_frag_del

    mut = np.random.choice(xrange(0, 4),
                           p=[p_frag_copy, p_frag_move, p_frag_del, p_nothing])

    offspring = bit_flip(muts[mut](genome), p_bit_flip)

    return offspring


def evaluate(ind, problem_db, lbl_db, nb_instances=-1):
    """Evaluate indiv NN on nb_instances
    with input=problem_db and output=lbl_db"""
    #TODO eval always on the same instances
    fit = 0
    if nb_instances == -1:
        instances = problem_db
        answers = lbl_db
    else:
        #Randomly sample "nbInstances" instances without replacement
        idx_instances = np.random.choice(xrange(len(problem_db)),
                                         nb_instances,
                                         replace=False)
        instances = [problem_db[instance] for instance in idx_instances]
        answers = [lbl_db[instance] for instance in idx_instances]

    for idx, inp in enumerate(instances):
        nn_answer = ind.activate(inp)
        #Error on each instance (euclidian distance to right answer)
        #is normalized between 0 and 1
        error = np.linalg.norm(np.subtract(nn_answer, answers[idx])) / \
            math.sqrt(len(inp))
        fit += 1.0 - error
    #Normalized for the size of the database:
    #average normalized euclidian error over instances
    return fit/float(len(instances))


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
    CONFIG = "config.ini"
    PARSER = ConfigParser.ConfigParser()
    PARSER.read(CONFIG)
    PARAMS = AutoVivification()
    VERBOSE = True

    for s in PARSER.sections():
        PARAMS[s] = map_config(s)
    print PARAMS
    sys.exit("")
    LOG_FILENAME = "logs/fitness_" + str(datetime.datetime.now()) + ".log"
    LOG_FILE = open(LOG_FILENAME, 'w')
###############################################################################
    N_OUT = 1  # One output for logical binary output
    N_IN = 2  # 8 for retina pb
    N_HID = 1
    NEUR_HID_LAYER = 2
    #Proportion of junk genes in-between genes on initialization of the genome
    FRAC_JUNK_GENES = 0.0
    SYMBOLS = ['0', '1']
    #Number of nucleotids for weight field on initialization
    #(the more nucleotids, the finer the resolution of initial weights)
    NB_WEIGHT_VALUES = 2
    NUCL_PER_WEIGHT = NB_WEIGHT_VALUES - 1
    #How to select start and end codons? Shorter (i.e. easier to randomly draw)
    #codons are more prone to disruptive variation? (are offspring viable?)
    #Grey-coding for links? minimizing "distances" between links?
    #Look for a symmetry in codons?
    #START_CODON = "111100111000110010010010"
    #END_CODON = "101011101101101100001111"
    #START_CODON = "111010010100"
    #END_CODON = "010100011011"
    START_CODON = "101001010011"
    END_CODON = "010110101100"
    #START_CODON = "11010"
    #END_CODON = "00110"
    PLEIOTROPY = True
    if N_HID == 0:
        N_LINKS = N_IN * N_OUT
    else:
        N_LINKS = N_IN * NEUR_HID_LAYER + (N_HID - 1) * NEUR_HID_LAYER**2 +\
            NEUR_HID_LAYER * N_OUT
    TGT_SIZE = int(math.ceil(max(1, math.log(N_LINKS, len(SYMBOLS)))))
#############################Evolutionary parameters###########################
    MU = 10
    NB_GENERATIONS = 15
    PARENT_CHILDREN_RATIO = 1.0  # number of children per parent
    # selectionRatio * mu == lambda nb of children
    LAMBDA = int(round(MU * PARENT_CHILDREN_RATIO))
    #mu parents in the beginning and lbda children every generation

###############################################################################
    PROBLEM_DB = []
    LABEL_DB = []
    PROBLEM_SIZE = N_IN
    DB_SIZE = -1  # -1 for whole problem
    TASK_SEQUENCE = ["t1.png", "t2.png"]
    # , "t3.png"]  # , "t4.png", "t5.png"]

    PROBLEM_ID = "IMG"
    # "RETAND" "MIN", "AND", "OR", "MAJ", "RETOR" "IMG"
###############################################################################
###############################################################################
##########################Evolutionary algorithm###############################
    INTERTASK_LOG = []
    POPULATION = [("", 0.0, [], "", None)] * MU

    #Initial genomes in population: mu individuals, valid ones
    #(all random links correctly encoded)
    for i in xrange(MU):
        g = create_genome(N_LINKS, FRAC_JUNK_GENES, NUCL_PER_WEIGHT,
                          (START_CODON, END_CODON), SYMBOLS)
        codons, PRETTY_GENOME_STRING = extract_codons(g, TGT_SIZE,
                                                      (START_CODON,
                                                       END_CODON),
                                                      pleio=PLEIOTROPY)
        net = map_to_standard_mlp(codons, N_IN, N_OUT, n_hidden=N_HID,
                                  neur_per_hid=NEUR_HID_LAYER)
        individual = (g, 0.0, codons, PRETTY_GENOME_STRING, net)
        POPULATION[i] = individual

###############################################################################
    #Total Time
    TIME_START = time.time()
    for each_task in TASK_SEQUENCE:
        POPULATION_LOG = []
        print "Learning task: ", each_task
        #Generate full problem
        if PROBLEM_ID == "IMG":
            IMG_FILENAME = each_task
            IMG = load_img(IMG_FILENAME)
            PROBLEM_DB = [[y / float(IMG[1]), x / float(IMG[0])]
                          for y in xrange(IMG[1]) for x in xrange(IMG[0])]
            LABEL_DB = [x / 255 for x in IMG[2]]
            TOTAL_INSTANCES = IMG[0] * IMG[1]
            # keep always same instances(todo)
        else:
            #Generate all instances
            TOTAL_INSTANCES = len(SYMBOLS)**PROBLEM_SIZE
            for i in xrange(TOTAL_INSTANCES):
                #binary in this case, different base for other symbol sets
                input_vector = [float(k) for k in
                                list(("{0:0%sb}" %
                                      str(PROBLEM_SIZE)).format(i))]
                PROBLEM_DB.append(input_vector)
                LABEL_DB.append(task(input_vector, each_task))

        INSTANCES_DB = []
        LABEL_INST_DB = []

        if DB_SIZE != -1:
            index_instances = np.random.choice(xrange(TOTAL_INSTANCES),
                                               DB_SIZE, replace=False)
            INSTANCES_DB.append([PROBLEM_DB[i] for i in index_instances])
            LABEL_INST_DB.append([LABEL_DB[i] for i in index_instances])
        else:
            INSTANCES_DB = PROBLEM_DB[:]
            LABEL_INST_DB = LABEL_DB[:]

        TIME_GEN_LOG = []

        PRETTY_GENOME_STRING = ""

        #Initial population: mu individuals, valid ones
        #(all random links correctly encoded) OR random binary string
        for i in xrange(MU):
            net = POPULATION[i][4]
            #Evaluate initial population
            fitness = evaluate(net, INSTANCES_DB, LABEL_INST_DB)
            POPULATION[i] = (POPULATION[i][0], fitness, POPULATION[i][2],
                             POPULATION[i][3], POPULATION[i][4])

        POPULATION_LOG.append(POPULATION)

        ITERATION = 0
        sys.stdout.write(str(ITERATION))

        #Evolutionary loop
        for ITERATION in xrange(1, NB_GENERATIONS + 1):
            #Time per generation
            time_gen_start = time.time()
            #Select lambda (possibly repeated) parents, rank-based
            parents = select(POPULATION, LAMBDA)

            children = [("", 0.0, [], "", None)] * LAMBDA

            #Generate lambda children by mutating lambda selected parents
            #Look for disruptive mutations (e.g. by measuring valid codons)
            #Note: it is easier to break a link gene than creating a new one
            for index, i in enumerate(parents):
                child = mutate(i[0][:])
                codons, PRETTY_GENOME_STRING = extract_codons(child[:],
                                                              TGT_SIZE,
                                                              (START_CODON,
                                                               END_CODON),
                                                              pleio=PLEIOTROPY)
                net = map_to_standard_mlp(codons, N_IN, N_OUT, n_hidden=N_HID,
                                          neur_per_hid=NEUR_HID_LAYER)
                #Evaluate children
                fitness = evaluate(net, INSTANCES_DB, LABEL_INST_DB)
                individual = (child[:], fitness, codons, PRETTY_GENOME_STRING,
                              net)
                children[index] = individual

            full_population = (POPULATION + children)
            #Truncation "plus" survivor sel./replacement [parents + children]
            surviving_individuals = survive(full_population, MU)

            POPULATION = [("", 0.0, [], "", None)] * MU
            for index, individual in enumerate(surviving_individuals):
                POPULATION[index] = individual

            POPULATION_LOG.append(POPULATION)

            #Logging operations at the end of generation
            sys.stdout.write(" - " + str(ITERATION))
            time_gen_end = time.time()
            TIME_GEN_LOG.append([time_gen_end - time_gen_start])

        INTERTASK_LOG.append(POPULATION_LOG)
        formatted_fitness = [[format_float(indiv[1]) for indiv in generation]
                             for generation in POPULATION_LOG]
        LOG_FILE.write(str(formatted_fitness) + "\n")
        #transfer to next task: best individual
        max_fitness = POPULATION[0][1]
        best_individual = POPULATION[0]
        for individual in POPULATION:
            if individual[1] > max_fitness:
                max_fitness = individual[1]
                best_individual = individual
        #If there is a following task, i.e. current is not last task
        if each_task != TASK_SEQUENCE[-1]:
            #Initialize population with mutated copies of best_individual
            for index in xrange(MU):
                altered_copy = mutate(best_individual[0][:])
                codons, PRETTY_GENOME_STRING = extract_codons(g, TGT_SIZE,
                                                              (START_CODON,
                                                               END_CODON),
                                                              pleio=PLEIOTROPY)
                net = map_to_standard_mlp(codons, N_IN, N_OUT, n_hidden=N_HID,
                                          neur_per_hid=NEUR_HID_LAYER)
                individual = (altered_copy, 0.0, codons,
                              PRETTY_GENOME_STRING, net)
                POPULATION[index] = individual

        print
###############################################################################
    LOG_FILE.close()
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
            codons, PRETTY_GENOME_STRING = extract_codons(g[0], TGT_SIZE,
                                                          (START_CODON,
                                                           END_CODON),
                                                          pleio=PLEIOTROPY)
            net = map_to_standard_mlp(codons, N_IN, N_OUT, n_hidden=N_HID,
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
