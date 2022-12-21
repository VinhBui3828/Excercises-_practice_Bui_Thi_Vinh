import re
import matplotlib.pyplot as plt
import numpy as np
import requests
from scipy import special as sp
import xml.etree.ElementTree as ET
import os


if __name__ == '__main__':
    # Скачивание файла с исходными данными
    link = 'https://jenyay.net/uploads/Student/Modelling/task_02.xml'
    filename = link.split('/')[-1] # Присваиваем имя файла из его ссылки
    web = requests.get(link) # Записывает в "web" данные из файла по ссылке
    file = open(filename, 'wb').write(web.content) # Создает и записывает файл

    # Вывод нужного варианта и присвоение значений
    in_file = web.text # Присвоение текста содержащегося в файле

    file = open(filename, 'r', encoding='utf-8') # Открытие файла для прочтения
    print (file.read()) # Чтение файла

    # Паттерн для поиска по тексту файла
    
    reg = r"variant number=\"3\" D=\"(\d+e-\d+|\d+\.\d+)\" fmin=\"(\d+.\d+e\d+|\d+e-\d+)\" fmax=\"(\d+e\d+|\d+\.\d+)\"/"

    matches = re.search(reg, in_file, re.MULTILINE) # Поиск по паттерну

    file.close() # Закрытие файла

    # ПРОВЕРКА !!!
    #print (matches.group(0))
    #print(D, fmin, fmax)

    # Присвоение значений из файла под нужным вариантом
    D = float(matches.group(1))
    fmin = float(matches.group(2))
    fmax = float(matches.group(3))
    
    r = D / 2
    c = 3*10**8

    # Функции для упрощения записи
    def hn(n, x):
        return (sp.jn(n, x) + 1j*sp.yn(n, x))

    def an(n, x):
        return (sp.jn(n, x) / hn(n, x))
    
    def bn(n, x):
        return ((k*r*sp.jn(n-1, x) - n*sp.jn(n,x)) /
                            (k*r*hn(n-1, x) - n*hn(n,x)))



    # Расчет значений
    f = np.linspace(fmin, fmax, 2500)
    lambda_ = c / f
    k = 2*np.pi / lambda_
    
    # Расчет суммы
    s = 0
    for n in range(1, 1000, 1):
        if all((-1)**n * (n + 0.5) * (bn(n, k*r) - an(n, k*r))) >= 10e-12:
            s += ((-1)**n * (n + 0.5) * (bn(n, k*r) - an(n, k*r)))
            #print(s)
        else:
            break

    sigma = ((lambda_**2/np.pi) * (np.abs(s))**2)

if __name__ == '__main__':
   
    xmin = 0.0
    xmax = 100.0
    count = 200
    plt.plot(f, sigma)
   
    plt.xlabel("$f$,Гц")
    plt.ylabel("$σ$,м^2")
    plt.grid()
    
    plt.show()



    




#!!!Для наглядности можно вывести в консоль
#print (f)
#print (lambda_)
#print(sigma)

    f = [str(j) for j in f]
    lambda_ = [str(j) for j in lambda_]
    sigma = [str(j) for j in sigma]

# Make file XML

if __name__ == '__main__':


   

    # Создание главного элемента:

    root = ET.Element('data')

    # Создание подэлементов:

    for i in range(len(sigma)):
        frequencydata = ET.SubElement(root, 'frequencydata')
        f = ET.SubElement(frequencydata, 'f').text = str(f)


        lambdadata = ET.SubElement(root, 'lambdadata')
        lambda_ = ET.SubElement(lambdadata, 'lambda_').text = str(lambda_)

        rcsdata = ET.SubElement(root, 'rcsdata')
        sigma = ET.SubElement(rcsdata, 'sigma').text = str(sigma)


    # Создание директории если её нет 
    d = os.path.dirname(__file__) # директория скрипта
    p = r'{}/results'.format(d) # создание path

    try:
        os.makedirs(p)
    except OSError:
        pass

    # Сохранение файла XML
    tree = ET.ElementTree(root)
    root  = tree.getroot()
    tree.write('results/results.xml', encoding='utf-8', xml_declaration=True)
    
   

