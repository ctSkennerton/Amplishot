from ctypes import *

class WuManber(object):
    def __init__(self,keys,so='_wumanber.so'):
        """ Initialise the WuManber object with required parameters
            Use __loadText__ and __loadKeywords__ to generate CTypes
            @keys:  list, string or filename
            @so:    name of the shared library linked to
        """
        self.so = CDLL(so)
        self.keywords =  None
        self.clist_of_cstrings = None # NOT A PYTHON TYPE
        self.len_clist_of_strings = None # NOT A PYTHON TYPE
        self.ctext = None # NOT A PYTHON TYPE
        self.len_ctext = None # NOT A PYTHON TYPE
        self.text = None
        self.wm = self.so.wumanber_init()
        self.__loadKeywords__(keys)
        
    def __loadText__(self,text):
        """ Parse the text provided by __init__. Depending on the type
            and whether the text is actually a URL, read the text into
            memory and create a CType c_char_p
            @text: string,url or filename
        """
        self.ctext = c_char_p(text)
        self.len_ctext = c_int(len(text))
    
      
    def __loadKeywords__(self,keys):
        """ Depending on the type() of keys, first create a Python list of
            keywords and then convert that to a C array of CType c_char_p
            @keys:  list, string or filename
        """
        self.keywords = keys
        self.clist_of_cstrings = (c_char_p*(len(self.keywords)))()
        self.len_clist_of_strings = c_int(len(self.clist_of_cstrings))
        i = 0
        for pystring in self.keywords:
            self.clist_of_cstrings[i] = c_char_p(pystring)
            i+=1
        self.so.wumanber_load_patterns(self.wm, self.clist_of_cstrings, self.len_clist_of_strings)
              
    def search_text(self,text):
        """ search_text is responsible for actually performing the text search
            @nocase:  boolean, whether to use case sensitive searching or not
            @verbose  boolean, whether to use the callback to print results or not
            @returns: int, the number of matches found in the text
        """

        self.__loadText__(text)
        found_pos,found_pattern_index = c_int(), c_int()
        self.so.wumanber_search_text(self.wm,self.len_ctext, self.ctext,self.clist_of_cstrings,
          byref(found_pos),byref(found_pattern_index))
        if found_pattern_index != -1:
            return found_pos.value, self.keywords[found_pattern_index.value]
        else:
            return None