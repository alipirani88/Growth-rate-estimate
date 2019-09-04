__author__ = 'alipirani'
import ConfigParser
Config = ConfigParser.ConfigParser()
#change the name of config file in case using different one.
Config.read("./config")


# Read Config File
def ConfigSectionMap(section, Config):
    dict1 = {}
    if not Config.has_section(section):
        #keep_logging('ERROR: Please Check the section name: \'{}\' in config file'.format(section), 'Please Check the section name: \'{}\' in config file'.format(section), logger, 'exception')
        print "ERROR: Please Check the section name: %s in config file" % section
        exit()
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1
