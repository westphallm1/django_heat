from django.shortcuts import render
from django.conf import settings
# Create your views here.
from django.http import HttpResponse

def home(request):
    return HttpResponse(settings.BASE_DIR)