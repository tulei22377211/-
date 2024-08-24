import cv2
import numpy as np
import matplotlib.pyplot as plt



img_clr_1 = cv2.imread('lena.png')
img_1 = cv2.cvtColor(img_clr_1, cv2.COLOR_BGR2GRAY)

img_clr_2 = cv2.imread('tian.jpg',cv2.IMREAD_COLOR)  # 读取彩色图像
img_clr_2 =  cv2.resize(img_clr_2, (int(img_clr_2.shape[1]*0.25),int(img_clr_2.shape[0]*0.25)),interpolation=cv2.INTER_CUBIC)
img_2 = cv2.cvtColor(img_clr_2, cv2.COLOR_BGR2GRAY)  # 转换为灰度图像



##对比度拉伸
def contrast_stretch(img):
    p2, p98 = np.percentile(img, (2, 98))  # 计算图像的2%和98%分位数
    img_rescale = cv2.convertScaleAbs(img, alpha=255.0/(p98-p2), beta=-p2*255.0/(p98-p2))  # 对比度拉伸
    img_rescale = img_rescale.astype(np.uint8)  # 转换数据类型
    return img_rescale

img_duibidu_1 = contrast_stretch(img_1)  # 对比度拉伸
img_duibidu_2 = contrast_stretch(img_2)  # 对比度拉伸


## 直方图均衡化
def hist_equalization(img):
    img_equalized = cv2.equalizeHist(img)  # 直方图均衡化
    return img_equalized

img_junhenghua_1 = hist_equalization(img_1)  # 直方图均衡化
img_junhenghua_2 = hist_equalization(img_2)  # 直方图均衡化


## 自适应直方图均衡化
def adapthist_equalization(img):
    clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8))  # 创建自适应直方图均衡化对象
    img_equalized = clahe.apply(img)  # 自适应直方图均衡化
    return img_equalized

img_zhishiyinjunhenghua_1 = adapthist_equalization(img_1)  # 自适应直方图均衡化
img_zhishiyinjunhenghua_2 = adapthist_equalization(img_2)  # 自适应直方图均衡化



images = [img_1,img_1,1, img_2,img_2,1,
          img_1, img_duibidu_1,1, img_2, img_duibidu_2,1,
          img_1, img_junhenghua_1,1, img_2, img_junhenghua_2,1,
          img_1, img_zhishiyinjunhenghua_1,1, img_2, img_zhishiyinjunhenghua_2,1,
          ]

titles = ['原图_1', '原图_1','直方图','原图_2', '原图_2','直方图',
          '原图_1', '对比度拉伸_1','直方图','原图_2', '对比度拉伸_2', '直方图',
          '原图_1', '均衡化_1','直方图','原图_2', '均衡化_2', '直方图',
          '原图_1', '自适应均衡化_1','直方图','原图_2', '自适应均衡化_2', '直方图',
          ]

plt.rcParams['font.sans-serif'] = ['SimHei']  # 或者其他CJK字体
plt.rcParams['axes.unicode_minus'] = False

x = 4   #x行
y = 6   #y列

for i in range(x):
    # 绘制原图_1
    plt.subplot(x, y, i * 6 + 1)
    plt.imshow(images[i * 6], 'gray', vmin=0, vmax=255)
    plt.title(titles[i * 6], fontsize=8)
    plt.xticks([]), plt.yticks([])

    #绘制结果
    plt.subplot(x, y, i * 6 + 2)
    plt.imshow(images[i * 6+1], 'gray', vmin=0, vmax=255)
    plt.title(titles[i * 6 + 1], fontsize=8)
    plt.xticks([]), plt.yticks([])

    # 绘制直方图 plt.hist,ravel 函数将数组降成一维
    plt.subplot(x, y, i * 6 + 3)
    plt.hist(images[i * 6+1].ravel(), 256, [0, 256])
    plt.title(titles[i * 6 + 2], fontsize=8)
    #plt.axvline(median, color='r', linestyle='--', linewidth=1)
    plt.xlim([0, 256])

    # 绘制原图_2
    plt.subplot(x, y, i * 6 + 4)
    plt.imshow(images[i * 6+3], 'gray', vmin=0, vmax=255)
    plt.title(titles[i * 6+3], fontsize=8)
    plt.xticks([]), plt.yticks([])

    #绘制结果
    plt.subplot(x, y, i * 6 + 5)
    plt.imshow(images[i * 6+4], 'gray', vmin=0, vmax=255)
    plt.title(titles[i * 6 + 4], fontsize=8)
    plt.xticks([]), plt.yticks([])

    # 绘制直方图 plt.hist,ravel 函数将数组降成一维
    plt.subplot(x, y, i * 6 + 6)
    plt.hist(images[i * 6+4].ravel(), 256, [0, 256])
    plt.title(titles[i * 6 + 5], fontsize=8)
    #plt.axvline(median, color='r', linestyle='--', linewidth=1)
    plt.xlim([0, 256])

    

    """ plt.subplot(2, 2, i * 2 + 2)
    for j,col in enumerate(color):
        histr = cv2.calcHist([images[i*2]],[j],None,[256],[0,256])
        plt.plot(histr,color = col)
        plt.xlim([0,256])
    plt.title(titles[i * 2 + 1], fontsize=8) """
plt.show()
