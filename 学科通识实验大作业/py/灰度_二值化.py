import cv2
import numpy as np
import matplotlib.pyplot as plt


img_clr_1 = cv2.imread('lena.png')
img_1 = cv2.cvtColor(img_clr_1, cv2.COLOR_BGR2GRAY)

img_clr_2 = cv2.imread('tian.jpg',cv2.IMREAD_COLOR)  # 读取彩色图像
img_clr_2 =  cv2.resize(img_clr_2, (int(img_clr_2.shape[1]*0.25),int(img_clr_2.shape[0]*0.25)),interpolation=cv2.INTER_CUBIC)
img_2 = cv2.cvtColor(img_clr_2, cv2.COLOR_BGR2GRAY)  # 转换为灰度图像


""" ## 非局部均值
def non_local_mean(img, h=3, templateWindowSize=7, searchWindowSize=21):
    img_filtered = cv2.fastNlMeansDenoising(img, None, h, templateWindowSize, searchWindowSize)  # 非局部均值
    return img_filtered

img_nonlocal = non_local_mean(img)  # 非局部均值 """



""" ## 均值平滑

def mean_smoothing(img, ksize=3):
    img_smoothed = cv2.blur(img, (ksize, ksize))  # 均值平滑
    return img_smoothed

img_mean = mean_smoothing(img)  # 均值平滑

# 显示图像
cv2.imshow('original', img)
cv2.imshow('smoothed', img_mean)
cv2.waitKey()
cv2.destroyAllWindows()  """


""" ## 高斯平滑

def gaussian_smoothing(img, ksize=3, sigma=0):
    img_smoothed = cv2.GaussianBlur(img, (ksize, ksize), sigma)  # 高斯平滑
    return img_smoothed

img_gaussian = gaussian_smoothing(img)  # 高斯平滑

# 显示图像
cv2.imshow('original', img)
cv2.imshow('smoothed', img_gaussian)
cv2.waitKey()
cv2.destroyAllWindows()  """


""" ## 中值滤波

def median_filtering(img, ksize=3):
    img_filtered = cv2.medianBlur(img, ksize)  # 中值滤波
    return img_filtered

img_median = median_filtering(img)  # 中值滤波

# 显示图像
cv2.imshow('original', img)
cv2.imshow('filtered', img_median)
cv2.waitKey()
cv2.destroyAllWindows() """


""" ## 双边滤波

def bilateral_filtering(img, d=9, sigmaColor=75, sigmaSpace=75):
    img_filtered = cv2.bilateralFilter(img, d, sigmaColor, sigmaSpace)  # 双边滤波
    return img_filtered

img_bilateral = bilateral_filtering(img)  # 双边滤波

# 显示图像
cv2.imshow('original', img)
cv2.imshow('filtered', img_bilateral)
cv2.waitKey()
cv2.destroyAllWindows() """



## 边缘检测
def edge_detection(img, ksize=3):
    blur = cv2.GaussianBlur(img, (ksize, ksize), 0)  # 高斯滤波
    img_edges = cv2.Canny(blur, 100, 200)  # 边缘检测
    contours, hierarchy = cv2.findContours(img_edges, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)  # 轮廓检测
    #cv2.drawContours(img, contours, -1, (0,255,0), 3)  # 绘制轮廓
    return img_edges

img_edges_1 = edge_detection(img_1)  # 边缘检测
img_edges_2 = edge_detection(img_2)  # 边缘检测

""" # 显示图像
cv2.imshow('original', img_1)
cv2.imshow('edges', img_edges_1)
cv2.imshow('original', img_2)
cv2.imshow('edges', img_edges_2)
cv2.waitKey()
cv2.destroyAllWindows() """


## sobel算子
def sobel_operator(img, ksize=3):
    #img_gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)  # 灰度化
    img_x = cv2.Sobel(img, cv2.CV_16S, 1, 0, ksize=ksize)  # x方向sobel算子
    img_y = cv2.Sobel(img, cv2.CV_16S, 0, 1, ksize=ksize)  # y方向sobel算子
    img_absx = cv2.convertScaleAbs(img_x)  # 取绝对值
    img_absy = cv2.convertScaleAbs(img_y)  # 取绝对值
    img_sobel = cv2.addWeighted(img_absx, 0.5, img_absy, 0.5, 0)  # 合并x,y方向的结果
    return img_sobel

img_sobel_1 = sobel_operator(img_1)  # sobel算子
img_sobel_2 = sobel_operator(img_2)  # sobel算子

""" # 显示图像
cv2.imshow('original', img_1)
cv2.imshow('sobel', img_sobel_1)
cv2.imshow('original', img_2)
cv2.imshow('sobel', img_sobel_2)
cv2.waitKey()
cv2.destroyAllWindows() """


## scharr算子
def scharr_operator(img):
    #img_gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)  # 灰度化
    img_x = cv2.Scharr(img, cv2.CV_16S, 1, 0)  # x方向scharr算子
    img_y = cv2.Scharr(img, cv2.CV_16S, 0, 1)  # y方向scharr算子
    img_absx = cv2.convertScaleAbs(img_x)  # 取绝对值
    img_absy = cv2.convertScaleAbs(img_y)  # 取绝对值
    img_scharr = cv2.addWeighted(img_absx, 0.5, img_absy, 0.5, 0)  # 合并x,y方向的结果
    return img_scharr

img_scharr_1 = scharr_operator(img_1)  # scharr算子
img_scharr_2 = scharr_operator(img_2)  # scharr算子

""" # 显示图像
cv2.imshow('original', img_1)
cv2.imshow('scharr', img_scharr_1)
cv2.imshow('original', img_2)
cv2.imshow('scharr', img_scharr_2)
cv2.waitKey()
cv2.destroyAllWindows() """


## prewitt算子
def prewitt_operator(img, ksize=3):
    #img_gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)  # 灰度化
    # 定义 Prewitt 算子的核
    kernel_x = np.array([[-1, 0, 1],
                        [-1, 0, 1],
                        [-1, 0, 1]])

    kernel_y = np.array([[-1, -1, -1],
                        [0, 0, 0],
                        [1, 1, 1]])

    # 将 Prewitt 算子的核应用于图像
    img_x = cv2.filter2D(img, -1, kernel_x)
    img_y = cv2.filter2D(img, -1, kernel_y)
    img_absx = cv2.convertScaleAbs(img_x)  # 取绝对值
    img_absy = cv2.convertScaleAbs(img_y)  # 取绝对值
    img_prewitt = cv2.addWeighted(img_absx, 0.5, img_absy, 0.5, 0)  # 合并x,y方向的结果
    return img_prewitt

img_prewitt_1 = prewitt_operator(img_1)  # prewitt算子
img_prewitt_2 = prewitt_operator(img_2)  # prewitt算子

""" # 显示图像
cv2.imshow('original', img_1)
cv2.imshow('prewitt', img_prewitt_1)
cv2.imshow('original', img_2)
cv2.imshow('prewitt', img_prewitt_2)
cv2.waitKey()
cv2.destroyAllWindows() """


## roberts算子
def roberts_operator(img, ksize=3):
    #img_gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)  # 灰度化
    def roberts_operator(img):
        roberts_x = np.array([[-1, 0], [0, 1]])
        roberts_y = np.array([[0, -1], [1, 0]])
        grad_x = cv2.filter2D(img, -1, roberts_x)
        grad_y = cv2.filter2D(img, -1, roberts_y)
        grad = np.sqrt(grad_x ** 2 + grad_y ** 2)
        return grad
    
    img_roberts = roberts_operator(img)  # roberts算子
    return img_roberts

img_roberts_1 = roberts_operator(img_1)  # roberts算子
img_roberts_2 = roberts_operator(img_2)  # roberts算子

img_roberts_1 = np.array(img_roberts_1, dtype=np.uint8)
img_roberts_2 = np.array(img_roberts_2, dtype=np.uint8)

""" # 显示图像
cv2.imshow('original', img_1)
cv2.imshow('roberts', img_roberts_1)
cv2.imshow('original', img_2)
cv2.imshow('roberts', img_roberts_2)
cv2.waitKey()
cv2.destroyAllWindows() """


## laplacian算子
def laplacian_operator(img, ksize=3):
    #img_gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)  # 灰度化
    img_laplacian = cv2.Laplacian(img, cv2.CV_16S, ksize=ksize)  # laplacian算子
    img_abs = cv2.convertScaleAbs(img_laplacian)  # 取绝对值
    return img_abs

img_laplacian_1 = laplacian_operator(img_1)  # laplacian算子
img_laplacian_2 = laplacian_operator(img_2)  # laplacian算子

""" # 显示图像
cv2.imshow('original', img_1)
cv2.imshow('laplacian', img_laplacian_1)
cv2.imshow('original', img_2)
cv2.imshow('laplacian', img_laplacian_2)
cv2.waitKey()
cv2.destroyAllWindows() """


## 固定阈值二值化
def thresholding(img):
    threshold=np.mean(img)
    img_binary = cv2.threshold(img, threshold, 255, cv2.THRESH_BINARY)[1]  # 固定阈值二值化
    return img_binary

img_binary_1 = thresholding(img_1)  # 固定阈值二值化
img_binary_2 = thresholding(img_2)  # 固定阈值二值化

""" # 显示图像
cv2.imshow('original', img_1)
cv2.imshow('binary', img_binary_1)
cv2.imshow('original', img_2)
cv2.imshow('binary', img_binary_2)
cv2.waitKey()
cv2.destroyAllWindows() """


## 自适应阈值二值化
def adaptive_thresholding(img, blockSize=11, C=2):
    img_binary = cv2.adaptiveThreshold(img, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, blockSize, C)  # 自适应阈值二值化
    return img_binary

img_adaptive_1 = adaptive_thresholding(img_1)  # 自适应阈值二值化
img_adaptive_2 = adaptive_thresholding(img_2)  # 自适应阈值二值化

""" # 显示图像
cv2.imshow('original', img_1)
cv2.imshow('adaptive', img_adaptive_1)
cv2.imshow('original', img_2)
cv2.imshow('adaptive', img_adaptive_2)
cv2.waitKey()
cv2.destroyAllWindows() """


## otsu二值化
def otsu_thresholding(img):
    ret, img_binary = cv2.threshold(img, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)  # otsu二值化
    return img_binary

img_otsu_1 = otsu_thresholding(img_1)  # otsu二值化
img_otsu_2 = otsu_thresholding(img_2)  # otsu二值化

""" # 显示图像
cv2.imshow('original', img_1)
cv2.imshow('otsu', img_otsu_1)
cv2.imshow('original', img_2)
cv2.imshow('otsu', img_otsu_2)
cv2.waitKey()
cv2.destroyAllWindows() """


## 局部otsu二值化
def local_otsu_thresholding(img, blockSize=11, C=2):
    img_blur = cv2.GaussianBlur(img, (5, 5), 0)  # 高斯滤波 
    cv2.imshow('blur', img_blur)   
    img_binary = cv2.adaptiveThreshold(img_blur, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, blockSize, C)  # 自适应阈值二值化
    return img_binary

img_localotsu_1 = local_otsu_thresholding(img_1)  # 局部otsu二值化
img_localotsu_2 = local_otsu_thresholding(img_2)  # 局部otsu二值化

""" # 显示图像
cv2.imshow('original', img_1)
cv2.imshow('localotsu', img_localotsu_1)
cv2.imshow('original', img_2)
cv2.imshow('localotsu', img_localotsu_2)
cv2.waitKey()
cv2.destroyAllWindows() """


images = [img_1,img_1,1, img_2,img_2,1,
          #img_1,img_edges_1,1, img_2,img_edges_2,1,
          #img_1,img_sobel_1,1, img_2,img_sobel_2,1,
          #img_1,img_laplacian_1,1, img_2,img_laplacian_2,1, 
          img_1,img_binary_1,1, img_2,img_binary_2,1,
          img_1,img_adaptive_1,1, img_2,img_adaptive_2,1,
          img_1,img_otsu_1,1, img_2,img_otsu_2,1,
          img_1,img_localotsu_1,1, img_2,img_localotsu_2,1,
          ]

titles = ['原图_1', '原图_1','直方图','原图_2', '原图_2','直方图',
          '原图_1', '固定阈值', '直方图', '原图_2', '固定阈值', '直方图',
          '原图_1', '自适应阈值', '直方图', '原图_2', '自适应阈值', '直方图',
          '原图_1', 'Otsu二值化', '直方图', '原图_2', 'Otsu二值化', '直方图',
          '原图_1', '局部Otsu二值化', '直方图', '原图_2', '局部Otsu二值化', '直方图',
          ]

plt.rcParams['font.sans-serif'] = ['SimHei']  # 或者其他CJK字体
plt.rcParams['axes.unicode_minus'] = False

x = 5   #x行
y = 4   #y列

for i in range(x):
    # 绘制原图_1
    plt.subplot(x, y, i * 4 + 1)
    plt.imshow(images[i * 6], 'gray', vmin=0, vmax=255)
    plt.title(titles[i * 6], fontsize=8)
    plt.xticks([]), plt.yticks([])

    #绘制结果
    plt.subplot(x, y, i * 4 + 2)
    plt.imshow(images[i * 6+1], 'gray', vmin=0, vmax=255)
    plt.title(titles[i * 6 + 1], fontsize=8)
    plt.xticks([]), plt.yticks([])

    """ # 绘制直方图 plt.hist,ravel 函数将数组降成一维
    plt.subplot(x, y, i * 6 + 3)
    plt.hist(images[i * 6+1].ravel(), 256, [0, 256])
    plt.title(titles[i * 6 + 2], fontsize=8)
    #plt.axvline(median, color='r', linestyle='--', linewidth=1)
    plt.xlim([0, 256]) """

    # 绘制原图_2
    plt.subplot(x, y, i * 4 + 3)
    plt.imshow(images[i * 6+3], 'gray', vmin=0, vmax=255)
    plt.title(titles[i * 6+3], fontsize=8)
    plt.xticks([]), plt.yticks([])

    #绘制结果
    plt.subplot(x, y, i * 4 + 4)
    plt.imshow(images[i * 6+4], 'gray', vmin=0, vmax=255)
    plt.title(titles[i * 6 + 4], fontsize=8)
    plt.xticks([]), plt.yticks([])

    """ # 绘制直方图 plt.hist,ravel 函数将数组降成一维
    plt.subplot(x, y, i * 6 + 6)
    plt.hist(images[i * 6+4].ravel(), 256, [0, 256])
    plt.title(titles[i * 6 + 5], fontsize=8)
    #plt.axvline(median, color='r', linestyle='--', linewidth=1)
    plt.xlim([0, 256]) """

    


plt.show()



















