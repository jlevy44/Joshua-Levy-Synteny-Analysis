import os
def generateConfigs(KaryotypeFiles,linkFile,configPath=os.getcwd()+'/',KaryotypePath=os.getcwd()+'/',LinkPath=os.getcwd()+'/',subgenomesCheck=0):
    """generateConfigs inputs
    configPath- path to export configuration files for Circos to
    KaryotypePath- location of karyotype files
    Link Path - location of link files
    LinkFile - filename of specific pairwise synteny comparison that a circos figure is being generated out of
    subgenomesCheck- if subgenome is being accessed, color set of links for chr a different color even though color for
    each ideogram/chromosome is the same
    NO return value for function, but generates configuration files circos.conf and linksAndrules.conf for particular
    pairwise comparison, to be used to generate circos plot"""
    # text for configuration file for circos.conf
    configText = """# circos.conf
    # 240:circos-0.69 jlevy$ /Applications/circos-0.69/bin/circos -conf /Applications/circos-0.69/bin/circos.conf

    karyotype = %s,%s

    chromosomes_units = 1000000
    #chromosomes_scale = /./=1rn
    chromosomes_display_default = yes

    # include ideogram configuration
    <<include %stxideogram.conf>>

    ## include ticks configuration
    <<include %stxticks.conf>>

    <colors>
    </colors>

    # include links and rules
    <<include %slinksAndrules.conf>>



    ################################################################
    # The remaining content is standard and required. It is imported
    # from default files in the Circos distribution.
    #
    # These should be present in every Circos configuration file and
    # overridden as required. To see the content of these files,
    # look in etc/ in the Circos distribution.

    <image>
    # Included from Circos distribution.
    <<include etc/image.conf>>
    </image>

    # RGB/HSV color definitions, color lists, location of fonts, fill patterns.
    # Included from Circos distribution.
    <<include etc/colors_fonts_patterns.conf>>

    # Debugging, I/O an dother system parameters
    # Included from Circos distribution.
    <<include etc/housekeeping.conf>>"""%(KaryotypePath+KaryotypeFiles[0],KaryotypePath+KaryotypeFiles[1],configPath,
                                          configPath,configPath)
    #^^^ above added in are the paths/filenames of the karyotype files, the configuration file paths in three different
    # spots

    # the links and rules configuration file text is split up into three parts
    linkAndruleChunk1 = """# links and rules

    <links>

    <link>
    file = %s
    radius = 0.99r
    bezier_radius = 0r
    #thickness = 2
    ribbon = yes
    color = black_a4

    <rules>
    <rule>
    condition = var(intrachr)
    show = no
    </rule>

    """%(LinkPath+linkFile) # link file and path added here to be accessed by circos

    # let's make  linkAndruleChunk2, import the rules for the link files, that is how to color the links
    # looks like...
    """<rule>
    condition = to(ChrSy)
    color = chr14
    </rule>"""
    # the coloring rules are generated from the colors on the karyotype files
    ruleInputFile = open(KaryotypePath+KaryotypeFiles[1],'r')
    ruleList = []
    if subgenomesCheck: # if subgenome like 524 or 523 being accessed, color set of links for chr a different color even
        # though color for each ideogram/chromosome is the same
        count = 1
        # create list of rules to be written to file
        for line in ruleInputFile:
            # string detailing how to color the links that come from certain target chromosomes for target species
            ruleList += ['\n<rule>\ncondition = to(%s)\ncolor = chr%d\n</rule>\n' % (line.split()[2],count)]
            count += 1
    else:
        # create list of rules to be written to file
        for line in ruleInputFile:
            # string detailing how to color the links that come from certain target chromosomes for target species
            lineList = line.split()
            ruleList += ['\n<rule>\ncondition = to(%s)\ncolor = %s\n</rule>\n'%(lineList[2],lineList[-1])]
    ruleInputFile.close()

    # compile the list of rules together
    linkAndruleChunk2 = ''.join(str(rule) for rule in ruleList)

    # third segment of link and rules chunk, ending syntax
    linkAndruleChunk3 = """

    </rules>
    </link>


    </links>"""

    # output for link and rule config file
    linkAndruleConfigOut = linkAndruleChunk1+linkAndruleChunk2+linkAndruleChunk3

    # generate Config Files circos and linksandrules
    open(configPath+'circos.conf','w').close()
    open(configPath+'linksAndrules.conf','w').close()

    # writing outputs to files
    circosConfig = open(configPath+'circos.conf','w')
    circosConfig.write(configText)
    circosConfig.close()

    linksAndRuleConfig = open(configPath+'linksAndrules.conf','w')
    linksAndRuleConfig.write(linkAndruleConfigOut)
    linksAndRuleConfig.close()