class Person(object):
    def __init__(self, firstname, lastname):
        self.firstname = firstname
        self.lastname = lastname
    def getFirstname(self):
        return self.firstname
    def getLastname(self):
        return self.lastname
    def getFullname(self):
        return f"{self.firstname} {self.lastname}"
        
class Student(Person):
    def __init__(self, firstname, lastname, subject):
        super().__init__(firstname, lastname)
        self.subject = subject
    
    def printNameSubject(self):
        print(f"{self.getFullname()}, {self.subject}")

