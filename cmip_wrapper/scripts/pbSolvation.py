#!/usr/bin/python3
#
# MDWeb Titration wrapper
import CMIP
def main(cen):
    gr=CMIP.Grid()
    gr.setreadGrid(2)
    gr.perfill=0.95
    gr.int = [0.5,0.5,0.5]

    gr0=CMIP.Grid(1);
    gr0.setreadGrid(2)
    gr0.perfill=0.7
    gr0.int[1.5,1.5,1.5]

    inpP=CMIP.InputParams("AA-Neutral Prot PB Solvation", gr, gr0)
    inpP.addKeyword(
    {   'tipcalc':0,
        'irest':0,
        'fullrst':0,
        'coorfmt':2,
        'calcgrid':1,
        'carmip':1,
        'orest':0,
        'rstonly':0,
        'pbelec':1,
        'pbinic':2,
        'pbfocus':1,
        'solvenergy':1
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
