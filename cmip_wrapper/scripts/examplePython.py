#!/usr/bin/python3

import CMIP

def main():
    gr = CMIP.Grid()
    gr.setreadGrid(2)
    gr.int = [0.5,0.5,0.5]
    gr0 = CMIP.Grid(1)
    gr0.setreadGrid(2)
    gr0.int = [1.5,1.5,1.5]
    inp = CMIP.InputParams("CMIP test",gr,gr0)
    inp.addKeyword({'CUBEOUTPUT':1,'MINOUT':1})
    run = CMIP.Run(inp)
    run.addFile({'hs': "1kim_l.pdb", 'cube' : "1kim_l.cube"})
    result = run.execute()
    print (result.asString())
    #print (result.getFileContents("cube"))


if __name__ == "__main__":
    main()
