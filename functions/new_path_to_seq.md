# New algo for path to seq. 

identify at each step which kind of folding step(nothing,simple addition, refolding, extension)

based on step kind add domains  

maybe for different steps different domain lengths? 

maybe keep track of all occuring substructures in form of a set --> check at each step if substructure allready happend

Right now 4 different kind of steps:



## Refolding

Simplest AFP with full refolding each step: 

```bash
.               |   b
()              |   b   b*
.()             |   b   b* a*   a  b
()()            |   b   b* a*   a  b c  c* b* a
.()()           |   b   b* a*   a  b c  c* b* a* d*  d a b c
()()()          |   b   b* a*   a  b c  c* b* a* d*  d a b c e  e* c* b* a* d* 
```




Addition of a module which does not get paired upon Transcription. Can only occur on uneven steps. 

- only few modifications 
- potential modifications at all module if they would have the same affinity to the new module  



## Extension
```bash
.      |     b   
()     |     b   b*             
.()    |     b   b* a*  a b     #refold event -> addition of domains
(())   |     b   b* a*  a b  b* #no mod necessary to achieve this path
```

Additional module which gets transcribed gets paired with module which was unpaired in the previous step. 

I think no modifications necessary

- Modifications not necessary since it binds with the open module anyways



## Nothing Step
```bash
.    |   b
()   |   b   b*
().  | a b   b* a  b   vs.   b a  a b*   b # unsure which variant better 
```


Addition of a module which does not get paired upon Transcription. Can only occur on uneven steps. 

- only few modifications 
- potential modifications at all module if they would have the same affinity to the new module  

