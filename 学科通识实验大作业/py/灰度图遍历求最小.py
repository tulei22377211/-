import cv2
import numpy as np
import matplotlib.pyplot as plt


def convert_value(x, a, b, c):
    if x < c:
        return a
    else:
        return b


def bianli_find_min_array(img,A,A_0):
    min_norm = float('inf')
    norm = float('inf')
    min_array = None
    min_a = 0
    min_b = 255

    #初值
    ret3, th3 = cv2.threshold(img, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    median = int(ret3)
    

    #min_a, min_b, median = 0, 255, np.median(A)
    #median = np.mean(A)
    
    min_norm = float('inf')
    step = 80
    print(min_a, min_b, median)

    
    for a in range(int(0), min(median,int(median+step)), 10):
        for b in range(max(median,int(median-step)), int(255), 10):
            converted_A = np.vectorize(convert_value)(A, a, b, median)
            norm = np.linalg.norm(converted_A - A_0)
            if norm < min_norm:
                min_norm = norm
                min_array = converted_A
                min_a, min_b = a, b
    for c in range(max(min_a,int(median-step)), min(min_b,int(median+step)), max(1,int(step/8))):
            converted_A = np.vectorize(convert_value)(A, min_a, min_b, c)
            norm = np.linalg.norm(converted_A - A_0)
            if norm < min_norm:
                min_norm = norm
                min_array = converted_A
                median = c
        
    print(min_a, min_b, median)
                
    while step > 0:
        for a in range(int(min_a-step), int(min_a+step), max(1,int(step/5))):
            for b in range(max(median,int(min_b-step)), int(min_b+step), max(1,int(step/5))):
                converted_A = np.vectorize(convert_value)(A, a, b, median)
                norm = np.linalg.norm(converted_A - A_0)
                if norm < min_norm:
                    min_norm = norm
                    min_array = converted_A
                    min_a, min_b = a, b
        for c in range(max(min_a,int(median-step)), min(min_b,int(median+step)), max(1,int(step/10))):
            converted_A = np.vectorize(convert_value)(A, min_a, min_b, c)
            norm = np.linalg.norm(converted_A - A_0)
            if norm < min_norm:
                min_norm = norm
                min_array = converted_A
                median = c
        print(min_a, min_b, median)

        # 如果当前norm和之前的norm相差小于1，则退出循环
        
        step //= 2

    for a in range(int(min_a-3), int(min_a+3), 1):
        for c in range(int(median-3), int(median+3), 1):
            for b in range(int(min_b-3), int(min_b+3), 1):
                converted_A = np.vectorize(convert_value)(A, a, b, c)
                norm = np.linalg.norm(converted_A - A_0)
                if norm < min_norm:
                    min_norm = norm
                    min_array = converted_A
                    min_a, min_b, median = a, b, c
    
    print(min_a, min_b, median)
                    
    return min_array,min_a,min_b,median,min_norm

""" min_b,norm = bianli_find_min_array(b)
min_g,norm = bianli_find_min_array(g)
min_r,norm = bianli_find_min_array(r) """

""" min_b = np.array(min_b, dtype=np.uint8)
min_g = np.array(min_g, dtype=np.uint8)
min_r = np.array(min_r, dtype=np.uint8)
min_array_rgb = cv2.merge((min_b,min_g,min_r)) """


# Read image
#img_clr = cv2.imread('lena.png')

img_clr = cv2.imread('tian.jpg')
img_clr =  cv2.resize(img_clr, (int(img_clr.shape[1]*0.25),int(img_clr.shape[0]*0.25)),interpolation=cv2.INTER_CUBIC)

# Convert to grayscale
img_gray = cv2.cvtColor(img_clr, cv2.COLOR_BGR2GRAY)

# 将图像转换为矩阵
A = np.array(img_gray, dtype=np.float32)

min_array,min_a,min_b,median,norm = bianli_find_min_array(img_gray,A,A)
min_array = np.array(min_array, dtype=np.uint8)

print('\n结束遍历')
print(norm)

cv2.imshow('Original', img_clr)
cv2.imshow('Grayscale', img_gray)
cv2.imshow('Min Array', min_array)
cv2.waitKey(0)

images = [img_gray, 0, min_array, 0]
titles = ['Original', 'Histogram', 
          'Min ', 'Histogram']

for i in range(2):
    # 绘制原图
    plt.subplot(2, 2, i * 2 + 1)
    plt.imshow(images[i * 2], 'gray', vmin=0, vmax=255)
    plt.title(titles[i * 2], fontsize=8)
    plt.xticks([]), plt.yticks([])

    # 绘制直方图 plt.hist,ravel 函数将数组降成一维
    plt.subplot(2, 2, i * 2 + 2)
    plt.hist(images[i * 2].ravel(), 256, [0, 256])
    plt.title(titles[i * 2 + 1], fontsize=8)
    plt.axvline(median, color='r', linestyle='--', linewidth=1)
    plt.xlim([0, 256])

    """ plt.subplot(2, 2, i * 2 + 2)
    for j,col in enumerate(color):
        histr = cv2.calcHist([images[i*2]],[j],None,[256],[0,256])
        plt.plot(histr,color = col)
        plt.xlim([0,256])
    plt.title(titles[i * 2 + 1], fontsize=8) """
plt.show()


""" cv2.imshow('Min Array RGB', min_array_rgb)
cv2.imshow('min_array r', min_r)
cv2.imshow('min_array g', min_g)
cv2.imshow('min_array b', min_b) """


