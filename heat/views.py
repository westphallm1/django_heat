from django.shortcuts import render
from django.conf import settings
from django.http import HttpResponse
from django.template import loader
import heat_cycle
import os

def home(request):
    #generate new image and solution table
    img_dir = os.path.join(settings.BASE_DIR,'static','default.png')
    txt_dir = os.path.join(settings.BASE_DIR,'static','default.txt')
    heat_cycle.gen_img_and_png(img_dir,txt_dir)
    
    #load the template
    template = loader.get_template('heat/home.html')
    img_url = os.path.join(settings.BASE_DIR,'static','default.png')
    context = {
        'root_dir': settings.BASE_DIR,
        'img_url': '/home/ubuntu/workspace/django_heat/static/default.png',
    }
    return HttpResponse(template.render(context, request))
    