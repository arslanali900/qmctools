import IOSection

class IOSectionClass:
    handle=0
    fileName=""
    def __init__(this):
        this.handle=IOSection.New()
    def OpenFile(this,fileName):
        return IOSection.OpenFile(this.handle,fileName)
    def GetName(this):
        return IOSection.GetName(this.handle)
    def NewFile(this,fileName):
        this.fileName=fileName
        return IOSection.NewFile(this.handle,fileName)
    
    def CloseFile(this):
        return IOSection.CloseFile(this.handle)
    def FlushFile(this):
        return IOSection.FlushFile(this.handle)
    def OpenSection2(this,name,num):
        return IOSection.OpenSectionNameNum(this.handle,name,num)
    def OpenSection(this,namenum):
        if type(namenum)==int:
            return IOSection.OpenSectionNum(this.handle,namenum)
        else:
            return IOSection.OpenSectionName(this.handle,namenum)
    def IncludeSection(this,name,fileName):
        return IOSection.IncludeSection(this.handle,name,fileName)
    def NewSection(this,name):
        return IOSection.NewSectionName(this.handle,name)
#        return IOSection.IncludeSection(this.handle,name,this.fileName)
    def CloseSection(this):
        IOSection.CloseSection(this.handle)
    def ReadVar(this,name):
        return IOSection.ReadVar(this.handle,name)
    def WriteVar(this,name,data):
        return IOSection.WriteVar(this.handle,name,data)
    def CountSections(this):
        return IOSection.CountSections(this.handle)
    def CountSections2(this,name):
        return IOSection.CountSectionsName(this.handle,name)
    def CountVars(this):
        return IOSection.CountVars(this.handle)
    def GetVarName(this,num):
        return IOSection.GetVarName(this.handle,num)
    def SetUnderscores(this,use):
        if (use):
            num = 1
        else:
            num = 0
        return IOSection.SetUnderscores(this.handle,num)


import os
class IOSectionClassList:
    def OpenFiles(this, baseName):
        this.IOlist = []
        proc=0
        done = 0
        while (done == 0):
            name = baseName + '.' + repr(proc) + '.h5'
            if (os.access(name,os.F_OK)):
                infile = IOSectionClass()
                success = infile.OpenFile (name)
                if (success==1) :
                    this.IOlist.append(infile)
            else:
                done=1
            proc = proc + 1
    def GetName(this):
        return this.IOlist[0].GetName()
    def OpenSection(this,namenum):
        success = 1
        for i in range(0,len(this.IOlist)):
            success = success and this.IOlist[i].OpenSection(namenum)
        return success
    def OpenSection2(this,name,num):
        success = 1
        for i in range(0,len(this.IOlist)):
            success = success and this.IOlist[i].OpenSection2(name,num)
        return success
    def CloseSection(this):
        map (lambda x: x.CloseSection(), this.IOlist)
    def ReadVar(this,name):
        data = map (lambda x: x.ReadVar(name),this.IOlist)
        return data
    def CountSections(this):
        return this.IOlist[0].CountSections()
    def CountSections2(this,name):
        return this.IOlist[0].CountSections2(name)
    def CountVars(this):
        return this.IOlist[0].CountVars()
    def GetVarName(this,num):
        return this.IOlist[0].GetVarName(num)
    def len(this):
        return len(this.IOlist)
    


    

#a=IOSectionClass()
#a.OpenFile("grr.h5")
#a.OpenSection("Energies")
#print len(a.ReadVar("PotentialEnergy"))
#a.CloseSection()
#a.OpenSection("ppPC");
#a.OpenSection("grid");
#r = a.ReadVar ("Points");
#print r

## b = IOSectionClass()
## b.OpenFile ("junk.dat")
## print b.ReadVar ("string_data")
## print b.ReadVar ("int_2");
## print b.ReadVar ("double_2");
## print b.ReadVar ("bool_2");
## print b.ReadVar ("string_2");

## print b.ReadVar ("int_3");
## print b.ReadVar ("double_3");
## print b.ReadVar ("bool_3");
## print b.ReadVar ("string_3");

## print 
## print b.ReadVar ("int_4");
## print 
## print b.ReadVar ("double_4");
## print 
## print b.ReadVar ("bool_4");
## print 
## print b.ReadVar ("string_4");
