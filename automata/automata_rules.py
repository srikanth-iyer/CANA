

## automata rules samples
## list of DCT automata
from cana.boolean_node import BooleanNode


automata_output_list = {
    'GKL'   : ['0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '1', '1', '1', '1'],
    'GP'    : ['0', '0', '0', '0', '0', '1', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '0', '1', '0', '1', '0', '0', '0', '0', '0', '1', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1'],
    'GEP_1' : ['0', '0', '0', '1', '0', '0', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '0', '0', '1', '0', '0', '0', '0', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '1', '1', '1', '1', '0', '0', '0', '1', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '0', '0', '1', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1'],
    'GEP_2' : ['0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '0', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '0', '1', '1', '1', '0', '0', '0', '0', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '1', '1', '1', '1', '0', '1', '1', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '1', '1', '0', '1', '1', '1'],
    'Davis' : ['0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '1', '1', '0', '1', '0', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '0', '0', '1', '1', '1', '1', '0', '0', '0', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '0', '0', '1', '0', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '0', '0', '1', '1', '1', '1', '1'],
    'Das'   : ['0', '0', '0', '0', '0', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '0', '0', '0', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '0', '0', '0', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '0', '0', '0', '1', '1', '1', '1', '0', '0', '1', '1', '0', '0', '0', '1', '0', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1'],
    'ABK'   : ['0', '0', '0', '0', '0', '1', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '0', '1', '0', '1', '0', '0', '0', '0', '0', '1', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1'],
    'DMC'   : ['0', '0', '0', '0', '0', '1', '0', '1', '0', '0', '0', '0', '0', '1', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '1', '0', '0', '0', '0', '1', '1', '1', '0', '0', '0', '0', '0', '1', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '1', '0', '1', '1', '1', '0', '1', '1', '1', '0', '0', '0', '0', '0', '0', '1', '1', '0', '1', '1', '1', '0', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '0', '0', '0', '0', '0', '1', '1', '0', '1', '1', '1', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '1', '0', '1', '1', '1', '0', '1', '1', '1', '1', '1', '1', '1'],
    'COE_1' : ['0', '0', '0', '0', '0', '0', '0', '1', '0', '0', '0', '1', '0', '1', '0', '0', '0', '0', '1', '1', '0', '0', '0', '0', '1', '1', '0', '1', '0', '1', '1', '1', '0', '0', '0', '1', '0', '0', '0', '1', '0', '0', '0', '0', '1', '1', '1', '1', '0', '0', '1', '1', '1', '0', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '0', '0', '0', '0', '0', '1', '0', '1', '1', '0', '1', '1', '0', '1', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '0', '0', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '0', '0', '0', '1', '0', '1', '1', '1', '1', '1', '0', '1', '1', '1', '1', '1', '1', '0', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1'],
    'COE_2' : ['0', '0', '0', '1', '0', '1', '0', '0', '0', '1', '0', '1', '0', '0', '0', '1', '0', '0', '1', '1', '0', '0', '0', '0', '0', '1', '0', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '0', '0', '0', '0', '1', '1', '0', '0', '1', '1', '1', '0', '0', '1', '0', '1', '1', '1', '1', '1', '0', '0', '0', '1', '0', '1', '1', '1', '0', '0', '0', '1', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '1', '1', '1', '1', '0', '0', '0', '0', '1', '1', '1', '1', '0', '1', '0', '1', '0', '0', '1', '1', '1', '1', '0', '0', '1', '1', '1', '1', '0', '1', '0', '1', '1', '1', '1', '1'],
    'MM401' : ['1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '0', '0', '1', '1', '1', '1', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '0', '0']
}


for output in automata_output_list:
    node = BooleanNode.from_output_list(outputs=automata_output_list[output])
    print(f"Bias of {output} = {node.bias()}")
## Effective graphs of the sample automata
# effective graphs of automata
automata = {
    "GKL": [["0###0#0", 0],["##0#0#0", 0],["###10#0", 0],["##00###", 0],["0##0###", 0],["1#1###1", 1],["###1##1", 1],["###11##", 1],["1#10###", 1],["1#1#1##", 1],],
    "GP": [["0#01###", 0],["0##10##", 0],["0#####0", 0],["###0##0", 0],["0#0#0##", 0],["##1#1#1", 1],["1#####1", 1],["1##1###", 1],["##10##1", 1],["###01#1", 1],],
    "GEP_1": [["000##0#", 0],["0###0#0", 0],["0##10##", 0],["0#0#00#", 0],["###0##0", 0],["00####0", 0],["##00#0#", 0],["00#1###", 0],["1#1###1", 1],["#11#1#1", 1],["#1#11##", 1],["1####11", 1],["#1##111", 1],["1##1###", 1],["##10##1", 1],["###0#11", 1],],
    "GEP_2": [["0#0###0", 0],["0####00", 0],["#0##000", 0],["###1#00", 0],["#00#0#0", 0],["##01##0", 0],["0##0###", 0],["#0#00##", 0],["11####1", 1],["111##1#", 1],["###1##1", 1],["11#0###", 1],["1##01##", 1],["1#1#11#", 1],["1###1#1", 1],["##11#1#", 1],],
    "Das": [["000##00", 0],["#010#00", 0],["0#01###", 0],["###00##", 0],["#001###", 0],["0#0#0##", 0],["##0#00#", 0],["##01#0#", 0],["00#0#00", 0],["#00#0##", 0],["##011#0", 0],["##1#11#", 1],["11#101#", 1],["###01#1", 1],["#1#01##", 1],["###011#", 1],["#11#1##", 1],["##1#1#1", 1],["11#1#11", 1],["1#001##", 1],["11##111", 1],["##11###", 1],],
    "Davis": [["#10#0#0", 0],["00#00##", 0],["##110#0", 0],["#00#0#1", 0],["10#011#", 0],["01##010", 0],["###1000", 0],["00#0#0#", 0],["00##000", 0],["#1#100#", 0],["001#0#0", 0],["0#1#010", 0],["##0#00#", 0],["#1#10#0", 0],["##00###", 0],["0##001#", 0],["#1#1#11", 1],["#011##1", 1],["0#1#11#", 1],["#110#0#", 1],["101#0#1", 1],["111##11", 1],["#001#10", 1],["101##01", 1],["1#1#10#", 1],["1110###", 1],["##11#11", 1],["#11#1##", 1],["1#100##", 1],["1#1#011", 1],["###11##", 1],["1#10#0#", 1],],
    "ABK": [["0#01###", 0],["0##10##", 0],["0#####0", 0],["###0##0", 0],["0#0#0##", 0],["##1#1#1", 1],["1#####1", 1],["1##1###", 1],["##10##1", 1],["###01#1", 1],],
    "DMC": [["00##01#", 0],["0#0###0", 0],["1011#01", 0],["01##000", 0],["#110001", 0],["1#00101", 0],["#01#010", 0],["#0#00#0", 0],["#111000", 0],["0#01#1#", 0],["1#10100", 0],["101110#", 0],["#0110#1", 0],["00##0#1", 0],["0##00##", 0],["##00000", 0],["#0##100", 0],["#010##0", 0],["00#0##0", 0],["#0000##", 0],["0#0#0##", 0],["0101###", 0],["#01101#", 0],["0##1100", 0],["#0#0#00", 0],["1000#0#", 0],["#00##00", 0],["01#1#00", 0],["###0111", 1],["1#0#11#", 1],["01101##", 1],["11###1#", 1],["110#0#1", 1],["#11#11#", 1],["#111#1#", 1],["0#1#1#1", 1],["1101###", 1],["1#01#1#", 1],["1##111#", 1],["##1#111", 1],["1#01##1", 1],["00##101", 1],["110#1#0", 1],["##101#1", 1],["1#10#11", 1],["#11#1#1", 1],["1###111", 1],["11#11##", 1],["1010##1", 1],["##1111#", 1],["0##01#1", 1],["#111##1", 1],["#001101", 1],["11100#0", 1],["11#1##1", 1],["#011000", 1],],
    "COE_1": [["01010##", 0],["00#01#0", 0],["#00#001", 0],["00#0#0#", 0],["#1#0110", 0],["01##000", 0],["#000##0", 0],["0#0#0#0", 0],["#00#1#0", 0],["0#00#0#", 0],["1#11#00", 0],["00101##", 0],["000###0", 0],["#0000##", 0],["##11010", 0],["0#00##0", 0],["#0#1100", 0],["#1#1000", 0],["0##0110", 0],["#00111#", 0],["00##100", 0],["##11100", 0],["1#0#110", 0],["#10010#", 0],["#1#0101", 0],["0##000#", 0],["1#110#0", 0],["101100#", 0],["0#0#00#", 0],["##001#0", 0],["0##1010", 0],["01#10#0", 0],["10#1001", 0],["#111#00", 0],["#1110#0", 0],["0##0#01", 0],["11###11", 1],["001100#", 1],["101#11#", 1],["110#01#", 1],["1##1011", 1],["01011##", 1],["10##101", 1],["##1#011", 1],["###1101", 1],["11##0#1", 1],["#1#0#11", 1],["#1#11#1", 1],["10#01#1", 1],["#11##11", 1],["##1111#", 1],["#111##1", 1],["101#1#1", 1],["01#111#", 1],["1#0101#", 1],["#1##111", 1],["1#100##", 1],["1#10#00", 1],["##111#1", 1],["10010#0", 1],["#0#1011", 1],["#110100", 1],["1010###", 1],["11#00##", 1],["##11#11", 1],["#10110#", 1],["11#1##1", 1],["1#1##11", 1],["##00111", 1],["1##0111", 1],["0#11##1", 1],["##1001#", 1],],
    "COE_2": [["00#011#", 0],["00#01#0", 0],["0#0###0", 0],["0##0111", 0],["001#11#", 0],["100#00#", 0],["##0110#", 0],["00#0#00", 0],["00101##", 0],["#1#001#", 0],["##01#00", 0],["00##110", 0],["01#0#11", 0],["1001#0#", 0],["#001##0", 0],["0100###", 0],["#1##010", 0],["##0#0#0", 0],["###10#0", 0],["00##000", 0],["00#000#", 0],["#1000##", 0],["0010#0#", 0],["0#0011#", 0],["##0000#", 0],["010#1##", 0],["00#1#10", 0],["#00##00", 0],["#11##01", 1],["1##1#11", 1],["##110#1", 1],["0##10#1", 1],["11#01##", 1],["1#10#0#", 1],["1#1#1##", 1],["#1#10#1", 1],["#1111##", 1],["#001#11", 1],["101###1", 1],["#11#10#", 1],["#11#1#0", 1],["#110#0#", 1],["10###11", 1],["11##11#", 1],["1010###", 1],["##1110#", 1],["1###111", 1],["#01001#", 1],["#0##011", 1],["#000101", 1],["##11#01", 1],["#111##1", 1],["###1011", 1],["1##01#1", 1],["1#11##1", 1],["1##011#", 1],["1#1##01", 1],],
    "MM401": [["11####1", 0],["1#1###1", 0],["###1##1", 0],["1#10###", 0],["##1111#", 0],["11#0###", 0],["1##01##", 0],["1#1#11#", 0],["1###1#1", 0],["0###0#0", 1],["0#0###0", 1],["#0000##", 1],["0####00", 1],["###10#0", 1],["###1#00", 1],["#00#0#0", 1],["##01##0", 1],["0##0###", 1],],
}
# for node in automata:
#     generated_node = BooleanNode.from_partial_lut(automata[node])
#     print(f"Bias of {node} = {generated_node.bias()}")

# Removing lines of the effective graph to make the automata incomplete, to regenerate the complete LUT from this and verify that the function works correctly. 


# removing some effective graph inputs to make them incomplete
incomplete_automata = {
    # "GKL": [["##0#0#0", 0],["###10#0", 0],["##00###", 0],["0##0###", 0],["###1##1", 1],["1#10###", 1],["1#1#1##", 1], ["1#1###1", 1],],  # ["0###0#0", 0], ["###11##", 1], is missing
    "GP": [["0#01###", 0],["0##10##", 0],["0#####0", 0],["###0##0", 0],["1#####1", 1],["1##1###", 1],["##10##1", 1],["###01#1", 1],], # ["0#0#00#", 0], is missing
    # "GEP_1": [["000##0#", 0],["0###0#0", 0],["0##10##", 0],["0#0#00#", 0],["###0##0", 0],["00####0", 0],["##00#0#", 0],["#11#1#1", 1],["#1#11##", 1],["1####11", 1],["#1##111", 1],["1##1###", 1],["##10##1", 1],["###0#11", 1],], # ["00#1###", 0],["1#1###1", 1], are missing
    # "GEP_2": [["0#0###0", 0],["0####00", 0],["#0##000", 0],["###1#00", 0],["#00#0#0", 0],["##01##0", 0],["0##0###", 0],["111##1#", 1],["###1##1", 1],["11#0###", 1],["1##01##", 1],["1#1#11#", 1],["1###1#1", 1],["##11#1#", 1],], # ["#0#00##", 0],["11####1", 1], are missing
    # "Das": [["000##00", 0],["#010#00", 0],["0#01###", 0],["###00##", 0],["#001###", 0],["0#0#0##", 0],["##0#00#", 0],["##01#0#", 0],["00#0#00", 0],["#00#0##", 0],["##011#0", 0],["##1#11#", 1],["11#101#", 1],["###01#1", 1],["#1#01##", 1],["###011#", 1],["#11#1##", 1],["##1#1#1", 1],["11#1#11", 1],["1#001##", 1],["11##111", 1],], # ["##11###", 1], is missing
    # "Davis": [["#10#0#0", 0],["00#00##", 0],["##110#0", 0],["#00#0#1", 0],["10#011#", 0],["01##010", 0],["###1000", 0],["00#0#0#", 0],["00##000", 0],["#1#100#", 0],["001#0#0", 0],["0#1#010", 0],["##0#00#", 0],["#1#10#0", 0],["##00###", 0],["0##001#", 0],["#1#1#11", 1],["#011##1", 1],["0#1#11#", 1],["#110#0#", 1],["101#0#1", 1],["111##11", 1],["#001#10", 1],["101##01", 1],["1#1#10#", 1],["1110###", 1],["##11#11", 1],["#11#1##", 1],["1#100##", 1],["1#1#011", 1],["###11##", 1],], # ["1#10#0#", 1], is missing
    "ABK": [["0#01###", 0],["0##10##", 0],["0#####0", 0],["###0##0", 0],["##1#1#1", 1],["1#####1", 1],["1##1###", 1],["##10##1", 1],], # ["0#0#0##", 0],["###01#1", 1], is missing
    # "DMC": [["00##01#", 0],["0#0###0", 0],["1011#01", 0],["01##000", 0],["#110001", 0],["#01#010", 0],["#0#00#0", 0],["#111000", 0],["0#01#1#", 0],["1#10100", 0],["101110#", 0],["#0110#1", 0],["00##0#1", 0],["0##00##", 0],["##00000", 0],["#0##100", 0],["#010##0", 0],["00#0##0", 0],["#0000##", 0],["0#0#0##", 0],["0101###", 0],["#01101#", 0],["0##1100", 0],["#0#0#00", 0],["1000#0#", 0],["#00##00", 0],["01#1#00", 0],["###0111", 1],["1#0#11#", 1],["01101##", 1],["11###1#", 1],["110#0#1", 1],["#11#11#", 1],["#111#1#", 1],["0#1#1#1", 1],["1101###", 1],["1#01#1#", 1],["1##111#", 1],["##1#111", 1],["1#01##1", 1],["00##101", 1],["110#1#0", 1],["##101#1", 1],["1#10#11", 1],["#11#1#1", 1],["1###111", 1],["11#11##", 1],["1010##1", 1],["##1111#", 1],["0##01#1", 1],["#111##1", 1],["#001101", 1],["11100#0", 1],["11#1##1", 1],], # ["1#00101", 0],["#011000", 1], is missing
    # "COE_1": [["01010##", 0],["00#01#0", 0],["#00#001", 0],["00#0#0#", 0],["#1#0110", 0],["01##000", 0],["#000##0", 0],["0#0#0#0", 0],["#00#1#0", 0],["0#00#0#", 0],["1#11#00", 0],["00101##", 0],["000###0", 0],["#0000##", 0],["##11010", 0],["0#00##0", 0],["#0#1100", 0],["#1#1000", 0],["0##0110", 0],["#00111#", 0],["00##100", 0],["##11100", 0],["1#0#110", 0],["#10010#", 0],["#1#0101", 0],["0##000#", 0],["1#110#0", 0],["101100#", 0],["0#0#00#", 0],["##001#0", 0],["0##1010", 0],["01#10#0", 0],["10#1001", 0],["#111#00", 0],["#1110#0", 0],["0##0#01", 0],["11###11", 1],["001100#", 1],["101#11#", 1],["110#01#", 1],["1##1011", 1],["01011##", 1],["10##101", 1],["##1#011", 1],["###1101", 1],["11##0#1", 1],["#1#0#11", 1],["#1#11#1", 1],["10#01#1", 1],["#11##11", 1],["##1111#", 1],["#111##1", 1],["101#1#1", 1],["01#111#", 1],["1#0101#", 1],["#1##111", 1],["1#100##", 1],["1#10#00", 1],["##111#1", 1],["10010#0", 1],["#0#1011", 1],["#110100", 1],["1010###", 1],["11#00##", 1],["##11#11", 1],["#10110#", 1],["11#1##1", 1],["1#1##11", 1],["##00111", 1],["1##0111", 1],["0#11##1", 1],], # ["##1001#", 1], is missing
    # "COE_2": [["00#011#", 0],["00#01#0", 0],["0#0###0", 0],["0##0111", 0],["001#11#", 0],["100#00#", 0],["##0110#", 0],["00#0#00", 0],["00101##", 0],["#1#001#", 0],["##01#00", 0],["00##110", 0],["01#0#11", 0],["1001#0#", 0],["#001##0", 0],["0100###", 0],["#1##010", 0],["##0#0#0", 0],["###10#0", 0],["00##000", 0],["00#000#", 0],["#1000##", 0],["0010#0#", 0],["0#0011#", 0],["##0000#", 0],["010#1##", 0],["00#1#10", 0],["#00##00", 0],["#11##01", 1],["1##1#11", 1],["##110#1", 1],["0##10#1", 1],["11#01##", 1],["1#10#0#", 1],["1#1#1##", 1],["#1#10#1", 1],["#1111##", 1],["#001#11", 1],["101###1", 1],["#11#10#", 1],["#11#1#0", 1],["#110#0#", 1],["10###11", 1],["11##11#", 1],["1010###", 1],["##1110#", 1],["1###111", 1],["#01001#", 1],["#0##011", 1],["#000101", 1],["##11#01", 1],["#111##1", 1],["###1011", 1],["1##01#1", 1],["1#11##1", 1],["1##011#", 1],], # ["1#1##01", 1], is missing
    # "MM401": [["11####1", 0],["1#1###1", 0],["###1##1", 0],["1#10###", 0],["##1111#", 0],["11#0###", 0],["1#1#11#", 0],["0#0###0", 1],["#0000##", 1],["0####00", 1],["###10#0", 1],["###1#00", 1],["#00#0#0", 1],["##01##0", 1], ["1###1#1", 0],["0##0###", 1], ], # ["1##01##", 0],["0###0#0", 1],are missing
}

automata_output_list = {
    'GKL'   : ['0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '1', '1', '1', '1'],
    'GP'    : ['0', '0', '0', '0', '0', '1', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '0', '1', '0', '1', '0', '0', '0', '0', '0', '1', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1'],
    # 'GEP_1' : ['0', '0', '0', '1', '0', '0', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '0', '0', '1', '0', '0', '0', '0', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '1', '1', '1', '1', '0', '0', '0', '1', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '0', '0', '1', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1'],
    # 'GEP_2' : ['0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '0', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '0', '1', '1', '1', '0', '0', '0', '0', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '1', '1', '1', '1', '0', '1', '1', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '1', '1', '0', '1', '1', '1'],
    # 'Davis' : ['0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '1', '1', '0', '1', '0', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '0', '0', '1', '1', '1', '1', '0', '0', '0', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '0', '0', '1', '0', '1', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '0', '0', '1', '1', '1', '1', '1'],
    # 'Das'   : ['0', '0', '0', '0', '0', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '0', '0', '0', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '0', '0', '0', '1', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '0', '0', '0', '1', '1', '1', '1', '0', '0', '1', '1', '0', '0', '0', '1', '0', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1'],
    'ABK'   : ['0', '0', '0', '0', '0', '1', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '0', '1', '0', '1', '0', '0', '0', '0', '0', '1', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1'],
    # 'DMC'   : ['0', '0', '0', '0', '0', '1', '0', '1', '0', '0', '0', '0', '0', '1', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '1', '0', '0', '0', '0', '1', '1', '1', '0', '0', '0', '0', '0', '1', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '1', '0', '1', '1', '1', '0', '1', '1', '1', '0', '0', '0', '0', '0', '0', '1', '1', '0', '1', '1', '1', '0', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '0', '0', '0', '0', '0', '1', '1', '0', '1', '1', '1', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '1', '0', '1', '1', '1', '0', '1', '1', '1', '1', '1', '1', '1'],
    # 'COE_1' : ['0', '0', '0', '0', '0', '0', '0', '1', '0', '0', '0', '1', '0', '1', '0', '0', '0', '0', '1', '1', '0', '0', '0', '0', '1', '1', '0', '1', '0', '1', '1', '1', '0', '0', '0', '1', '0', '0', '0', '1', '0', '0', '0', '0', '1', '1', '1', '1', '0', '0', '1', '1', '1', '0', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '0', '0', '0', '0', '0', '1', '0', '1', '1', '0', '1', '1', '0', '1', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '0', '0', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '0', '0', '0', '1', '0', '1', '1', '1', '1', '1', '0', '1', '1', '1', '1', '1', '1', '0', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1'],
    # 'COE_2' : ['0', '0', '0', '1', '0', '1', '0', '0', '0', '1', '0', '1', '0', '0', '0', '1', '0', '0', '1', '1', '0', '0', '0', '0', '0', '1', '0', '1', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '0', '0', '0', '0', '1', '1', '0', '0', '1', '1', '1', '0', '0', '1', '0', '1', '1', '1', '1', '1', '0', '0', '0', '1', '0', '1', '1', '1', '0', '0', '0', '1', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '1', '1', '1', '1', '0', '0', '0', '0', '1', '1', '1', '1', '0', '1', '0', '1', '0', '0', '1', '1', '1', '1', '0', '0', '1', '1', '1', '1', '0', '1', '0', '1', '1', '1', '1', '1'],
    # 'MM401' : ['1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '0', '1', '0', '0', '0', '1', '1', '1', '1', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '1', '0', '1', '0', '0', '0']
}

for automata in incomplete_automata:
    node = None
    generated_node_permuations = None
    node = BooleanNode.from_partial_lut(incomplete_automata[automata])
    '''print(f"Count of '?' {node.outputs.count('?')}")
    print(f"Count of '1' {node.outputs.count('1')}")
    print(f"Count of '0' {node.outputs.count('0')}")
    # indices where the output is '?'
    indices = [i for i, x in enumerate(node.outputs) if x == "?"]
    print(f"Indices of '?' {indices}")'''

    generated_node_permuations = BooleanNode.generate_with_required_bias( node, required_node_bias=0.5, verbose=True)
    list_of_output_lists = [node.outputs for node in generated_node_permuations]

    # print(automata)

    #checking if original rule is contained in the list of generated rules, to verify if the function is working correctly
    if automata_output_list[automata] in list_of_output_lists:
        print("Found a match")
    else:
        print("No match found")


    # # find the list in list_of_output_lists that matches automata_output_list[automata]
    # print("Index : ", list_of_output_lists.index(automata_output_list[automata]))