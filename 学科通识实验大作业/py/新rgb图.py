import cv2
import numpy as np
import matplotlib.pyplot as plt


## sobel算子
def sobel_operator(img, ksize=3):
    #img_gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)  # 灰度化
    img_x = cv2.Sobel(img, cv2.CV_16S, 1, 0, ksize=ksize)  # x方向sobel算子
    img_y = cv2.Sobel(img, cv2.CV_16S, 0, 1, ksize=ksize)  # y方向sobel算子
    img_absx = cv2.convertScaleAbs(img_x)  # 取绝对值
    img_absy = cv2.convertScaleAbs(img_y)  # 取绝对值
    img_sobel = cv2.addWeighted(img_absx, 0.5, img_absy, 0.5, 0)  # 合并x,y方向的结果
    return cv2.bitwise_not(img_sobel)  # 取反


## scharr算子
def scharr_operator(img):
    #img_gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)  # 灰度化
    img_x = cv2.Scharr(img, cv2.CV_16S, 1, 0)  # x方向scharr算子
    img_y = cv2.Scharr(img, cv2.CV_16S, 0, 1)  # y方向scharr算子
    img_absx = cv2.convertScaleAbs(img_x)  # 取绝对值
    img_absy = cv2.convertScaleAbs(img_y)  # 取绝对值
    img_scharr = cv2.addWeighted(img_absx, 0.5, img_absy, 0.5, 0)  # 合并x,y方向的结果
    return img_scharr


## laplacian算子
def laplacian_operator(img, ksize=3):
    #img_gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)  # 灰度化
    img_laplacian = cv2.Laplacian(img, cv2.CV_16S, ksize=ksize)  # laplacian算子
    img_abs = cv2.convertScaleAbs(img_laplacian)  # 取绝对值
    return cv2.bitwise_not(img_abs)  # 取反


##对比度拉伸
def contrast_stretch(img):
    p2, p98 = np.percentile(img, (2, 98))  # 计算图像的2%和98%分位数
    img_rescale = cv2.convertScaleAbs(img, alpha=255.0/(p98-p2), beta=-p2*255.0/(p98-p2))  # 对比度拉伸
    img_rescale = img_rescale.astype(np.uint8)  # 转换数据类型
    return img_rescale


## 自适应直方图均衡化
def adapthist_equalization(img):
    clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8))  # 创建自适应直方图均衡化对象
    img_equalized = clahe.apply(img)  # 自适应直方图均衡化
    return img_equalized


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
        for a in range(int(min_a-step), min(median,int(median+step)), max(1,int(step/5))):
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


def transfer(img):
    img_con_stretch = contrast_stretch(img)  # 对比度拉伸
    img_adp_equalized = adapthist_equalization(img)  # 自适应直方图均衡化

    img_sobel_adp = sobel_operator(img_adp_equalized)  # sobel算子
    img_scharr_adp = cv2.bitwise_not(scharr_operator(img_adp_equalized))  # scharr算子
    img_laplacian_adp = laplacian_operator(img_adp_equalized)  # laplacian算子

    img_sobel_con = sobel_operator(img_con_stretch)  # sobel算子
    img_scharr_con = cv2.bitwise_not(scharr_operator(img_con_stretch))  # scharr算子
    img_laplacian_con = laplacian_operator(img_con_stretch)  # laplacian算子

    
    img_lunkuo = cv2.addWeighted(img_sobel_adp, 0.5, img_scharr_adp, 0.5, 0)  # 融合sobel和laplacian算子
    img_jieguo = cv2.addWeighted(img_lunkuo, 0.2, img_adp_equalized, 0.8, 0)  # 融合融合后的图像和对比度拉伸后的图像

    # 显示结果
    #cv2.imshow('Original', img)
    #cv2.imshow('lunkuo', img_lunkuo)
    #cv2.imshow('laplacian_adp', img_laplacian_adp)
    #cv2.imshow('sobel_adp', img_sobel_adp)
    #cv2.imshow('scharr_adp', img_scharr_adp)

    #cv2.imshow('laplacian_con', img_laplacian_con)
    #cv2.imshow('sobel_con', img_sobel_con)
    #cv2.imshow('scharr_con', img_scharr_con)
    #cv2.imshow('jieguo', img_jieguo)
    cv2.waitKey(0)

    return img_jieguo


# Read image

#img_clr = cv2.imread('a.jpg')

#img_clr = cv2.imread('b.jpg')
#img_clr =  cv2.resize(img_clr, (int(img_clr.shape[1]*0.25),int(img_clr.shape[0]*0.25)),interpolation=cv2.INTER_CUBIC)

#img_clr = cv2.imread('d.png')
#img_clr =  cv2.resize(img_clr, (int(img_clr.shape[1]*0.25),int(img_clr.shape[0]*0.25)),interpolation=cv2.INTER_CUBIC)

img_clr = cv2.imread('lena.png')

#img_clr = cv2.imread('tian.jpg')
#img_clr =  cv2.resize(img_clr, (int(img_clr.shape[1]*0.25),int(img_clr.shape[0]*0.25)),interpolation=cv2.INTER_CUBIC)

#通道分割
b_0,g_0,r_0=cv2.split(img_clr)    #b = img[:, :, 0]

b = transfer(b_0)
g = transfer(g_0)
r = transfer(r_0)


print("\n开始遍历b_0")
min_b,min_b_a,min_b_b,min_b_c,min_b_norm = bianli_find_min_array(b_0,b_0,b_0)
print("\n开始遍历g_0")
min_g,min_g_a,min_g_b,min_g_c,min_g_norm = bianli_find_min_array(g_0,g_0,g_0)
print("\n开始遍历r_0")
min_r,min_r_a,min_r_b,min_r_c,min_r_norm = bianli_find_min_array(r_0,r_0,r_0)

print('\n开始遍历b')
min_b_trans,min_b_a_trans,min_b_b_trans,min_b_c_trans,min_b_norm_trans = bianli_find_min_array(b,b,b_0)
print('\n开始遍历g')
min_g_trans,min_g_a_trans,min_g_b_trans,min_g_c_trans,min_g_norm_trans = bianli_find_min_array(g,g,g_0)
print('\n开始遍历r')
min_r_trans,min_r_a_trans,min_r_b_trans,min_r_c_trans,min_r_norm_trans = bianli_find_min_array(r,r,r_0)

min_b = np.array(min_b, dtype=np.uint8)
min_g = np.array(min_g, dtype=np.uint8)
min_r = np.array(min_r, dtype=np.uint8)

min_b_trans = np.array(min_b_trans, dtype=np.uint8)
min_g_trans = np.array(min_g_trans, dtype=np.uint8)
min_r_trans = np.array(min_r_trans, dtype=np.uint8)

min_array_rgb = cv2.merge((min_b,min_g,min_r))
min_array_rgb_trans = cv2.merge((min_b_trans,min_g_trans,min_r_trans))

print('\n结束遍历')


#min_array,norm = bianli_find_min_array(A)
#min_array = np.array(min_array, dtype=np.uint8)

#print(norm)

cv2.imshow('Original', img_clr)
#cv2.imshow('Min Array', min_array)
cv2.imshow('Min Array RGB', min_array_rgb)
cv2.imshow('Min Array RGB Trans', min_array_rgb_trans)
""" cv2.imshow('min_array r', min_r)
cv2.imshow('min_array g', min_g)
cv2.imshow('min_array b', min_b) """
cv2.waitKey(0)

#cv2.imwrite('a_min_rgb.jpg', min_array_rgb)
#cv2.imwrite('a_trans_rgb.jpg', min_array_rgb_trans)

#cv2.imwrite('b_min_rgb.jpg', min_array_rgb)
#cv2.imwrite('b_trans_rgb.jpg', min_array_rgb_trans)

#cv2.imwrite('d_min_rgb.jpg', min_array_rgb)
#cv2.imwrite('d_trans_rgb.jpg', min_array_rgb_trans)

cv2.imwrite('lena_min_rgb.jpg', min_array_rgb)
cv2.imwrite('lena_trans_rgb.jpg', min_array_rgb_trans)

#cv2.imwrite('tian_min_rgb.jpg', min_array_rgb)
#cv2.imwrite('tian_trans_rgb.jpg', min_array_rgb_trans)

color = ['b','g','r']

images = [img_clr[:,:,[2,1,0]], 0, min_array_rgb[:,:,[2,1,0]], 0, min_array_rgb_trans[:,:,[2,1,0]], 0]

# 绘制图像
titles = ['Original', 'Histogram', 
          'Min RGB', 'Histogram',
          'Min RGB Trans', 'Histogram']

medians = [min_b_c, min_g_c, min_r_c]

for i in range(3):
    # 绘制原图
    plt.subplot(3, 2, i * 2 + 1)
    plt.imshow(images[i * 2], 'gray')
    plt.title(titles[i * 2], fontsize=8)
    plt.xticks([]), plt.yticks([])

    """ # 绘制直方图 plt.hist,ravel 函数将数组降成一维
    plt.subplot(3, 3, i * 3 + 2)
    plt.hist(images[i * 3].ravel(), 256)
    plt.title(titles[i * 3 + 1], fontsize=8)
    plt.xticks([]), plt.yticks([]) """

    plt.subplot(3, 2, i * 2 + 2)
    for j,col in enumerate(reversed(color)):
        histr = cv2.calcHist([images[i*2]],[j],None,[256],[0,256])
        plt.plot(histr,color = col)
        plt.axvline(medians[j],color = color[j],linestyle='--')
        plt.xlim([0,256])
    plt.title(titles[i * 2 + 1], fontsize=8)

plt.show()



""" for i,col in enumerate(color):
    histr = cv2.calcHist([img_clr],[i],None,[256],[0,256])
    plt.plot(histr,color = col)
    plt.xlim([0,256])
plt.show()

for i,col in enumerate(color):
    histr = cv2.calcHist([min_array_rgb],[i],None,[256],[0,256])
    plt.plot(histr,color = col)
    plt.xlim([0,256])

plt.show() """


