from functions import ir_domain_translator 
from functions import convert_functions as cv

d_seq_1 = "b l1 b* a* l2 a b c d l3 d* c* b*"

path_1 = [[".."],["(.).."],["..((.))..."],["(.)...(((.)))"]]

#extendet_path = cv.extended_fp_path(path,d_seq)



d_seq_2 = "b l b* a* l a b c d l d* c* b* l a b e f l f* e* b*" #bbbbblllllBBBBBAAADDDlllllaaabbbbbccclllllCCCBBBBBllllldddaaabbbbbeeelllllEEEBBBBB

d_seq_3 = "b l b* a* d* l a b c l c* b* l d a b e l e* b* " 
path_2 = [[".."],["(.)..."],["..((..)).."],["(.)....((.))."],["..(((..((.)).))).."],["(.)....((.))...((.))"]]


cv.extended_domain_path([d_seq_3])

d_seq_4 = "b l b* a* l a b l b* c* l c b d e l e* d* b*"

path_4 = [[".."],["(.).."],["..((.))."],["(.((.)).).."],["..((.)).((.))..."],["(.((.)).)...(((.)))"]]


d_seq_5 = "b a l b* c* l b d l a* b* l c b l d* b*"
path_5 = [["..."],["(..).."],["(..)....."],["((.(..)..))."],["...((.(...).))."],["(.)..((..(..).))"]]


d_seq_6 = "b d l b* a* l a b c l c* b* a* e* l e a b c l d* b*"

path_6 = [["..."],["(..).."],["...((.)).."],["(..)..(((.))).."],["(..)..(((.)))......."],["((.((.))..((((.)))).))"]]




ir_domain_translator.rna_design(seq= d_seq_6, path=path_6, out= 5)


 