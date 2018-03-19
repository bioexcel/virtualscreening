#!/usr/bin/python3
#
# MDWeb Titration wrapper
import CMIP
def main(cen):
    gr=CMIP.Grid()
    gr.setreadGrid(0);
    gr.int = [0.5,0.5,0.5]
    gr.cen = cen
    gr.dim = [40,40,40]

    gr0=CMIP.Grid(1)
    gr0.setreadGrid(2);
    gr0.int = [1.5,1.5,1.5]
    gr0.perfill= 0.6

    inpP=CMIP.InputParams("AA-Neutral Prot Interaction", gr, gr0)
    inpP.addKeyword(
    {   'tipcalc':3,
        'irest':0,
        'orest':0,
        'cutvdw':8.,
        'cutelec':10.,
        'coorfmt:>2,
        'calcgrid':1,
        'dields':2,
        'outegrid':1,
        'fvdw':0.8,    
        'pbelec':1
    })
    
    return inpP
        

#    gr = CMIP.Grid()
#    gr.setreadGrid(2)
#    gr.int = [0.5,0.5,0.5]
#    gr0 = CMIP.Grid(1)
#    gr0.setreadGrid(2)
#    gr0.int = [1.5,1.5,1.5]
#    inp = CMIP.InputParams("CMIP test",gr,gr0)
#    inp.addKeyword({'CUBEOUTPUT':1,'MINOUT':1})
#    run = CMIP.Run(inp)
#    run.addFile({'hs': "1kim_l.pdb", 'cube' : "1kim_l.cube"})
#    result = run.execute()
#    print (result.asString())
#    #print (result.getFileContents("cube"))


if __name__ == "__main__":
    main(argv)
