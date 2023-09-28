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


path_7 = [['..'], ['(.)..'], ['..((.))..'], ['(.)..(((.)))..'], ['..((.))..((((.))))..'], ['(.)..(((.)))..(((((.)))))..'], ['..((.))..((((.))))..((((((.))))))..'], ['(.)..(((.)))..(((((.)))))..(((((((.)))))))..'], ['..((.))..((((.))))..((((((.))))))..((((((((.))))))))..'], ['(.)..(((.)))..(((((.)))))..(((((((.)))))))..(((((((((.)))))))))..'], ['..((.))..((((.))))..((((((.))))))..((((((((.))))))))..((((((((((.)))))))))).'], ['(.)..(((.)))..(((((.)))))..(((((((.)))))))..(((((((((.))))))))).......(.....)']]

d_seq_7 = "b   l   b* a*  l  a b c  l  c* b* a* d*  l  d a b c e  l  e* c* b* a* d* f*  l  f d a b c   e g  l  g* e* c* b* a* d* f* h*  l  h f d a b c   e   g i  l  i* g* e* c* b* a* d* f* h* j*  l  j h f d a b c   e   g   i  l   b* "


path_8_steps = [['..'], ['(.)..'], ['..((.))..'], ['(.)..(((.)))..'], ['..((.))..((((.))))..'], ['(.)..(((.)))..(((((.)))))..'], ['..((.))..((((.))))..((((((.))))))..'], ['(.)..(((.)))..(((((.)))))..(((((((.)))))))..']]

d_seq_8_steps = "b   l   b* a*  l  a b c  l  c* b* a* d*  l  d a b c e  l  e* c* b* a* d* f*  l  f d a b c   e g  l  g* e* c* b* a* d* f* h*  l  h f d a b c   e   g i  l "
    

d_seq_6_steps = "b   l   b* a*  l  a b c  l  c* b* a* d*  l  d a b c e  l  e* c* b* a* d* f*"
path_6_steps = [['..'], ['(.)..'], ['..((.))..'], ['(.)..(((.)))..'], ['..((.))..((((.))))..'], ['(.)..(((.)))..(((((.))))).']]



"""
For testing nothing step 

ex: 
[["."],["()"],["()."],["(())"]]

seq1: a b c b* a* b c* b* a* 
seq2: b a c a* b* b c* a* b* 
"""

d_seq_8 = "a b c l b* a* l b l c* b* a*" #seq1 
d_path_8 = [["...."],["((..))."],["((..))..."],["(((.(..).)))"]]


d_seq_9 = "b a c l a* b* l b l c* a* b*" #seq2  
d_path_9 = [["...."],["((..))."],["((..))..."],["(((..(.).)))"]]


d_seq_10 = "c b a l a* b* l b l a* b* c*"
d_path_10 = [["...."],[".((.))."],[".((.))..."],["(((..(.).)))"]]

ir_domain_translator.rna_design(seq= d_seq_6_steps, path=path_6_steps, out= 6)


 