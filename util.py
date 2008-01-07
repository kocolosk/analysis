def getRun(path):
   """searches for an integer runnumber in the supplied path"""
   import re
   regex = re.search('(6|7)\d{6}',path)
   return int(regex.group())
