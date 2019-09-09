# ============================= #
# TESTING PYTHON CLASS - HOW TO #
# ============================= #

class Empty() :
  pass

class Name() :
  name = '"name" is the name of the attribute for the class Name()'

class Person() :
  def __init__(self, name, age):
    self.name = name
    self.age = age

# person = Person('Ambra', 26)
# print(person.name)
# print(person.age)
